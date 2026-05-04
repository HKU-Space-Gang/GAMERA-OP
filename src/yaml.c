#include "yaml.h"

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printError(const char *format, ...) {
  va_list args;
  va_start(args, format);
  fprintf(stderr, "\033[1;31m");
  vfprintf(stderr, format, args);
  fprintf(stderr, "\033[0m\n");
  va_end(args);
}

// Initialize a YAML document
void yaml_document_init(YamlDocument *doc) {
  doc->entries = NULL;
  doc->count = 0;
  doc->capacity = 0;
}

// Free a YAML document
void yaml_document_free(YamlDocument *doc) {
  for (size_t i = 0; i < doc->count; i++) {
    free(doc->entries[i].key);
    if (doc->entries[i].type == YAML_TYPE_PAIR) {
      free(doc->entries[i].value);
    } else if (doc->entries[i].type == YAML_TYPE_LIST) {
      for (size_t j = 0; j < doc->entries[i].list.count; j++) {
        free(doc->entries[i].list.items[j].value);
      }
      free(doc->entries[i].list.items);
    }
  }
  free(doc->entries);
  doc->entries = NULL;
  doc->count = 0;
  doc->capacity = 0;
}

// Add a key-value pair to the document
void yaml_document_add_pair(YamlDocument *doc, const char *key,
                            const char *value) {
  if (doc->count >= doc->capacity) {
    doc->capacity = doc->capacity == 0 ? 1 : doc->capacity * 2;
    doc->entries = realloc(doc->entries, doc->capacity * sizeof(YamlEntry));
  }
  doc->entries[doc->count].key = strdup(key);
  doc->entries[doc->count].type = YAML_TYPE_PAIR;
  doc->entries[doc->count].value = strdup(value);
  doc->count++;
}

// Add a list to the document
void yaml_document_add_list(YamlDocument *doc, const char *key) {
  if (doc->count >= doc->capacity) {
    doc->capacity = doc->capacity == 0 ? 1 : doc->capacity * 2;
    doc->entries = realloc(doc->entries, doc->capacity * sizeof(YamlEntry));
  }
  doc->entries[doc->count].key = strdup(key);
  doc->entries[doc->count].type = YAML_TYPE_LIST;
  doc->entries[doc->count].list.items = NULL;
  doc->entries[doc->count].list.count = 0;
  doc->entries[doc->count].list.capacity = 0;
  doc->count++;
}

// Add an item to the last list in the document
void yaml_document_add_list_item(YamlDocument *doc, const char *value) {
  if (doc->count == 0 || doc->entries[doc->count - 1].type != YAML_TYPE_LIST) {
    printError("Error: No list to add item to.");
    return;
  }

  YamlEntry *last_entry = &doc->entries[doc->count - 1];
  if (last_entry->list.count >= last_entry->list.capacity) {
    last_entry->list.capacity =
        last_entry->list.capacity == 0 ? 1 : last_entry->list.capacity * 2;
    last_entry->list.items =
        realloc(last_entry->list.items,
                last_entry->list.capacity * sizeof(YamlListItem));
  }
  last_entry->list.items[last_entry->list.count].value = strdup(value);
  last_entry->list.count++;
}

// Parse a YAML file into a document
int yaml_parse_file(YamlDocument *doc, const char *filename) {
  FILE *file = fopen(filename, "r");
  if (!file) {
    printError("Failed to open file: %s", filename);
    return -1;
  }

  char line[256];
  while (fgets(line, sizeof(line), file)) {
    // Trim leading whitespace
    char *ptr = line;
    while (*ptr == ' ') ptr++;

    // Skip empty or commented lines
    if (*ptr == '\n' || *ptr == '\0' || *ptr == '#' || *ptr == '\r' ||
        *ptr == '\t' || *ptr == '\f' || *ptr == '\v')
      continue;

    // Check if the line is a list item
    if (*ptr == '-') {
      ptr++;                      // Skip '-'
      while (*ptr == ' ') ptr++;  // Skip spaces after '-'
      char *value = strtok(ptr, "\n");
      if (value) {
        yaml_document_add_list_item(doc, value);
      }
    } else {
      // Parse key-value pair or list key
      char *key = strtok(ptr, ":");
      char *value = strtok(NULL, "\n");
      if (key && value) {
        while (*value == ' ') value++;  // Trim leading spaces in value
        yaml_document_add_pair(doc, key, value);
      } else if (key) {
        // If there's no value, assume it's the start of a list
        yaml_document_add_list(doc, key);
      }
    }
  }

  fclose(file);
  return 0;
}

// Get a string value by key
const char *yaml_get_string(YamlDocument *doc, const char *key) {
  for (size_t i = 0; i < doc->count; i++) {
    if (doc->entries[i].type == YAML_TYPE_PAIR &&
        strcmp(doc->entries[i].key, key) == 0) {
      return doc->entries[i].value;
    }
  }
  return NULL;  // Key not found
}

// Get an integer value by key
int yaml_get_int(YamlDocument *doc, const char *key, int *out_value) {
  const char *value = yaml_get_string(doc, key);
  if (!value) return 1;  // Key not found

  char *endptr;
  errno = 0;
  long result = strtol(value, &endptr, 10);
  if (errno != 0 || (*endptr != '\0' && *endptr != '\n' && *endptr != '\r'))
    return -1;  // Conversion failed
  *out_value = (int)result;
  return 0;
}

// Get a double value by key
//int yaml_get_double(YamlDocument *doc, const char *key, double *out_value) {
//  const char *value = yaml_get_string(doc, key);
//  if (!value) return 1;  // Key not found
//
//  char *endptr;
//  errno = 0;
//  double result = strtod(value, &endptr);
//  printf("key=%s, value=%s, result=%f\n", key, value, result);
//  if (errno != 0 || (*endptr != '\0' && *endptr != '\n' && *endptr != '\r'))
//    return -1;  // Conversion failed
//
//  *out_value = result;
//  return 0;
//}

int yaml_get_double(YamlDocument *doc, const char *key, double *out_value) {
  const char *value = yaml_get_string(doc, key);
  if (!value) return 1;

  char *endptr;
  double result = strtod(value, &endptr);

  // 跳过所有空白字符（空格、制表符、换行等）
  while (*endptr == ' ' || *endptr == '\t' || *endptr == '\n' || *endptr == '\r') {
    endptr++;
  }

  // 检查是否到达字符串结尾
  if (*endptr != '\0') {
    return -1;
  }

  *out_value = result;
  return 0;
}

// Get a list of strings by key
// int yaml_get_list_of_strings(YamlDocument *doc, const char *key,
//                              char ***out_list, size_t *out_count) {
//   for (size_t i = 0; i < doc->count; i++) {
//     if (doc->entries[i].type == YAML_TYPE_LIST &&
//         strcmp(doc->entries[i].key, key) == 0) {
//       *out_count = doc->entries[i].list.count;
//       *out_list = malloc(*out_count * sizeof(char *));
//       if (!*out_list) return -1;  // Memory allocation failed

//       for (size_t j = 0; j < *out_count; j++) {
//         (*out_list)[j] = strdup(doc->entries[i].list.items[j].value);
//         if (!(*out_list)[j]) return -1;  // Memory allocation failed
//       }
//       return 0;  // Success
//     }
//   }
//   return 1;  // Key not found
// }
int yaml_get_strings(YamlDocument *doc, const char *key, char **out_list,
                     size_t count) {
  for (size_t i = 0; i < doc->count; i++) {
    if (doc->entries[i].type == YAML_TYPE_LIST &&
        strcmp(doc->entries[i].key, key) == 0) {
      if (doc->entries[i].list.count != count) return -1;  // List size mismatch

      for (size_t j = 0; j < count; j++) {
        out_list[j] = strdup(doc->entries[i].list.items[j].value);
        if (!out_list[j]) return -1;  // Memory allocation failed
      }
      return 0;  // Success
    }
  }
  return 1;  // Key not found
}

// Get a list of integers by key
// int yaml_get_list_of_ints(YamlDocument *doc, const char *key, int **out_list,
// size_t *out_count) {
//     for (size_t i = 0; i < doc->count; i++) {
//         if (doc->entries[i].type == YAML_TYPE_LIST &&
//         strcmp(doc->entries[i].key, key) == 0) {
//             *out_count = doc->entries[i].list.count;
//             *out_list = malloc(*out_count * sizeof(int));
//             if (!*out_list) return -1; // Memory allocation failed

//             for (size_t j = 0; j < *out_count; j++) {
//                 char *endptr;
//                 errno = 0;
//                 long result = strtol(doc->entries[i].list.items[j].value,
//                 &endptr, 10); if (errno != 0 || *endptr != '\0') return -1;
//                 // Conversion failed
//                 (*out_list)[j] = (int)result;
//             }
//             return 0; // Success
//         }
//     }
//     return 1; // Key not found
// }

int yaml_get_ints(YamlDocument *doc, const char *key, int *out_list,
                  size_t count) {
  for (size_t i = 0; i < doc->count; i++) {
    if (doc->entries[i].type == YAML_TYPE_LIST &&
        strcmp(doc->entries[i].key, key) == 0) {
      if (doc->entries[i].list.count != count) return -1;  // List size mismatch
      for (size_t j = 0; j < count; j++) {
        char *endptr;
        errno = 0;
        long result = strtol(doc->entries[i].list.items[j].value, &endptr, 10);
        if (errno != 0 ||
            (*endptr != '\0' && *endptr != '\n' && *endptr != '\r'))
          return -1;  // Conversion failed
        out_list[j] = (int)result;
      }
      return 0;  // Success
    }
  }
  return 1;  // Key not found
}

// Get a list of doubles by key
// int yaml_get_list_of_doubles(YamlDocument *doc, const char *key, double
// **out_list, size_t *out_count) {
//     for (size_t i = 0; i < doc->count; i++) {
//         if (doc->entries[i].type == YAML_TYPE_LIST &&
//         strcmp(doc->entries[i].key, key) == 0) {
//             *out_count = doc->entries[i].list.count;
//             *out_list = malloc(*out_count * sizeof(double));
//             if (!*out_list) return -1; // Memory allocation failed

//             for (size_t j = 0; j < *out_count; j++) {
//                 char *endptr;
//                 errno = 0;
//                 double result = strtod(doc->entries[i].list.items[j].value,
//                 &endptr); if (errno != 0 || *endptr != '\0') return -1; //
//                 Conversion failed
//                 (*out_list)[j] = result;
//             }
//             return 0; // Success
//         }
//     }
//     return 1; // Key not found
// }

int yaml_get_doubles(YamlDocument *doc, const char *key, double *out_list,
                     size_t count) {
  for (size_t i = 0; i < doc->count; i++) {
    if (doc->entries[i].type == YAML_TYPE_LIST &&
        strcmp(doc->entries[i].key, key) == 0) {
      if (doc->entries[i].list.count != count) return -1;  // List size mismatch
      for (size_t j = 0; j < count; j++) {
        char *endptr;
        errno = 0;
        double result = strtod(doc->entries[i].list.items[j].value, &endptr);
        if (errno != 0 ||
            (*endptr != '\0' && *endptr != '\n' && *endptr != '\r'))
          return -1;  // Conversion failed
        out_list[j] = result;
      }
      return 0;  // Success
    }
  }
  return 1;  // Key not found
}
