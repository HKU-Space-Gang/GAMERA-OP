#ifndef GAMERA_NONORTHOGONAL_ABSOLUTE_CADENCE_H
#define GAMERA_NONORTHOGONAL_ABSOLUTE_CADENCE_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Target/source-independent coupling clock.  Due instants are integer
 * multiples of interval_s from physical t=0.  It contains no RC state and
 * performs no physics update; callers poll, perform a transactional refresh,
 * and commit only after that refresh succeeds.
 */
typedef struct gamera_no_absolute_cadence {
  double interval_s;
  double time_norm_s;
  double next_due_s;
} gamera_no_absolute_cadence;

/*
 * Initialize at a code-unit time.  require_alignment is intended for a
 * coupling whose held state is recomputed rather than checkpointed: a
 * mid-interval restart is rejected, while a restart at a due instant is
 * deterministic.
 */
int gamera_no_absolute_cadence_init(
    double interval_s, double time_norm_s, double initial_time_code,
    int require_alignment, gamera_no_absolute_cadence *cadence);

/* Clip a positive code-unit step so it lands exactly on next_due_s. */
int gamera_no_absolute_cadence_limit_timestep(
    const gamera_no_absolute_cadence *cadence, double time_code,
    double proposed_dt_code, double *limited_dt_code);

/*
 * Canonicalize a completed step before any cadence poll/source commit.
 * Snap a roundoff-close landing to the current absolute due instant, then to
 * an exactly coincident integration stop.  This is called before poll/source
 * commit; gaps larger than roundoff are unchanged.
 */
int gamera_no_absolute_cadence_canonicalize_step_end(
    const gamera_no_absolute_cadence *cadence, double time_code,
    double stop_time_code, double *canonical_time_code);

/*
 * Inspect a completed step without changing the clock.  Crossing a due
 * instant is an error; landing on it returns due=1 and its canonical physical
 * time.  This two-phase API keeps the clock unchanged if a source refresh
 * fails.
 */
int gamera_no_absolute_cadence_poll(
    const gamera_no_absolute_cadence *cadence, double time_code, int *due,
    double *due_time_s);

/* Advance by one interval after a successful refresh at due_time_s. */
int gamera_no_absolute_cadence_commit(
    gamera_no_absolute_cadence *cadence, double due_time_s);

/* Full scheduler restart tuple; restore validates without partial update. */
int gamera_no_absolute_cadence_export(
    const gamera_no_absolute_cadence *cadence, double *interval_s,
    double *time_norm_s, double *next_due_s);
int gamera_no_absolute_cadence_restore(
    double interval_s, double time_norm_s, double next_due_s,
    double restart_time_code, gamera_no_absolute_cadence *cadence);

int gamera_no_absolute_cadence_time_is_aligned(
    double time_code, double time_norm_s, double interval_s, int *aligned);

#ifdef __cplusplus
}
#endif

#endif
