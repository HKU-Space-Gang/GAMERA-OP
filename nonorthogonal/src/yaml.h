// yaml_parser.h
#ifndef YAML_PARSER_H
#define YAML_PARSER_H

#include <stddef.h>

void printError(const char *format, ...);

typedef enum { YAML_TYPE_PAIR, YAML_TYPE_LIST } YamlEntryType;

typedef struct {
  char *value;
} YamlListItem;

typedef struct {
  char *key;
  YamlEntryType type;
  union {
    char *value;  // For key-value pairs
    struct {
      YamlListItem *items;
      size_t count;
      size_t capacity;
    } list;  // For lists
  };
} YamlEntry;

typedef struct {
  YamlEntry *entries;
  size_t count;
  size_t capacity;
} YamlDocument;

// Initialize a YAML document.
void yaml_document_init(YamlDocument *doc);

// Free a YAML document and its associated memory.
void yaml_document_free(YamlDocument *doc);

// Add a key-value pair to the document.
void yaml_document_add_pair(YamlDocument *doc, const char *key,
                            const char *value);

// Add a list to the document.
void yaml_document_add_list(YamlDocument *doc, const char *key);

// Add an item to the last list in the document.
void yaml_document_add_list_item(YamlDocument *doc, const char *value);

// Parse a YAML file into a document.
int yaml_parse_file(YamlDocument *doc, const char *filename);

// Get a string value by key.
const char *yaml_get_string(YamlDocument *doc, const char *key);

// Get an integer value by key.
int yaml_get_int(YamlDocument *doc, const char *key, int *out_value);

// Get a double value by key.
int yaml_get_double(YamlDocument *doc, const char *key, double *out_value);

// Get a list of strings by key.
int yaml_get_strings(YamlDocument *doc, const char *key, char **out_list,
                     size_t count);

// Get a list of integers by key.
int yaml_get_ints(YamlDocument *doc, const char *key, int *out_list,
                  size_t count);

// Get a list of doubles by key.
int yaml_get_doubles(YamlDocument *doc, const char *key, double *out_list,
                     size_t count);

#endif  // YAML_PARSER_H
