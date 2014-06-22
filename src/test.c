#include <gep.h>
#include <visualize.h>
#include <utils.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <pthread.h>
#include <GL/glut.h>
#include <GL/freeglut_ext.h>

#define MAX_FRAMERATE 2

#define STRINGIZE(a) #a

#define CHECK_FIRST_LINE_INT(token) \
  if (0 == strncmp(STRINGIZE(token), line, strlen(STRINGIZE(token)))) \
  { \
    sscanf(line, "%s %d", token_buf, &short_params->token); \
  }

#define CHECK_NEXT_LINE(token) \
  else if (0 == strncmp(STRINGIZE(token), line, strlen(STRINGIZE(token))))

#define CHECK_NEXT_LINE_INT(token) \
  else if (0 == strncmp(STRINGIZE(token), line, strlen(STRINGIZE(token)))) \
  { \
    sscanf(line, "%s %d", token_buf, &short_params->token); \
  }

#define CHECK_NEXT_LINE_DBL(token) \
  else if (0 == strncmp(STRINGIZE(token), line, strlen(STRINGIZE(token)))) \
  { \
    sscanf(line, "%s %lf", token_buf, &short_params->token); \
  }

#define CHECK_NEXT_LINE_STR(token, container) \
  else if (0 == strncmp(STRINGIZE(token), line, strlen(STRINGIZE(token)))) \
  { \
    sscanf(line, "%s %s", token_buf, container); \
  }



static gep_ctx_t          *gep_ctx             = NULL;
static gep_stat_params_t   stat_params;

static visualize_ctx_t    *gr_ctx              = NULL;

static visualize_plot_t   *gr_plot_mean        = NULL;
static visualize_plot_t   *gr_plot_best        = NULL;
static visualize_plot_t   *gr_plot_target      = NULL;
static visualize_plot_t   *gr_plot_gep         = NULL;
static visualize_label_t  *gr_lb_progress      = NULL;

static unsigned long       simulation_started  = 0;
static unsigned long       last_frame_rendered = 0;
static unsigned long       last_output_printed = 0;

static volatile int        graph_enabled       = 0;
static volatile int        graph_to_stop       = 0;
static volatile int        graph_stopped       = 0;
static int                 progress            = 0;


/*----------------------------------------------------------------------------*/
//!
static void gep_sigterm_handler(int sig)
{
  if (SIGUSR1 == sig)
  {
    GEP_Stop(gep_ctx);
  }
}


/*----------------------------------------------------------------------------*/
//!
static void ui_KeyboardHandler(unsigned char key, int x, int y)
{
  switch (key)
  {
  case 27: //ESC
    GEP_Stop(gep_ctx);
    break;

  default:
    break;
  }

  /* Anti-Warning */
  x = x;
  y = y;
}


/*----------------------------------------------------------------------------*/
//!
static void stat_callback(double best_fitness, double avg_fitness, double best_mse, int generation, double *target_points, double *predicted_points)
{
  if (getTimeDiffMs(last_output_printed, getTimeNs()) > 5)
  {
    int ms_passed = (int)getTimeDiffMs(simulation_started, getTimeNs());
    printf("| %d\t\t%d\t\t%5.3lf\n", ms_passed, generation, best_mse);
    last_output_printed = getTimeNs();
  }

  progress = (100 * generation) / stat_params.total_generations_count;

  if (graph_enabled)
  {
    VISUALIZE_PlotAdd2dPoint(gr_plot_best, best_fitness);
    VISUALIZE_PlotAdd2dPoint(gr_plot_mean, avg_fitness);

    if (NULL != predicted_points)
    {
      VISUALIZE_PlotSetPoints(gr_plot_gep, predicted_points, stat_params.input_samples_count);
    }

    if (NULL != target_points)
    {
      VISUALIZE_PlotSetPoints(gr_plot_target, target_points, stat_params.input_samples_count);
    }
  }
}


/*----------------------------------------------------------------------------*/
//!
static void graph_redraw(void)
{
  if (graph_to_stop)
  {
    VISUALIZE_DestroyContext(gr_ctx);
    glutLeaveMainLoop();
    graph_stopped = 1;
    return;
  }

  if (graph_enabled && getTimeDiffMs(last_frame_rendered, getTimeNs()) > (1000 / MAX_FRAMERATE))
  {
    char progress_str[256];

    sprintf(progress_str, "Progress: %d%%", progress);
    VISUALIZE_LabelSetText(gr_lb_progress, progress_str);
    VISUALIZE_Render(gr_ctx);
    last_frame_rendered = getTimeNs();
  }
}


/*----------------------------------------------------------------------------*/
//!
static void graph_init(int argc, char *argv[])
{
  glutInit(&argc, argv);
  glutIdleFunc(graph_redraw);

  gr_ctx = VISUALIZE_CreateContext();

  glutKeyboardFunc(ui_KeyboardHandler);

  double gr_elements_basic_x_offset = 100;
  double gr_elements_basic_y_offset = 30;

  double label_x         = gr_elements_basic_x_offset;
  double label_y         = gr_elements_basic_y_offset;

  double plot_fit_width  = 200;
  double plot_fit_height = 200;

  double plot_best_x     = gr_elements_basic_x_offset;
  double plot_best_y     = label_y + 70;

  double plot_mean_x     = plot_best_x + plot_fit_width + 100;
  double plot_mean_y     = plot_best_y;



  double plot_f_width    = 200;
  double plot_f_height   = 200;
  double func_max        = 2.0;

  double plot_gep_x      = gr_elements_basic_x_offset;
  double plot_gep_y      = plot_best_y + plot_fit_height + 150;

  double plot_target_x   = plot_gep_x;
  double plot_target_y   = plot_gep_y + plot_f_height + 150;


  int points_per_dimension =
      (1 == stat_params.dimensions_count)
    ? (stat_params.input_samples_count)
    : ((int)sqrt(stat_params.input_samples_count));

  gr_plot_target = VISUALIZE_AddPlot(
    gr_ctx,
    plot_target_x, plot_target_y,
    plot_f_width,  plot_f_height,
    func_max,
    points_per_dimension,
    stat_params.dimensions_count + 1,
    "Target");

  gr_plot_gep = VISUALIZE_AddPlot(
    gr_ctx,
    plot_gep_x,   plot_gep_y,
    plot_f_width, plot_f_height,
    func_max,
    points_per_dimension,
    stat_params.dimensions_count + 1,
    "Predicted");

  gr_plot_best = VISUALIZE_AddPlot(
    gr_ctx,
    plot_best_x,    plot_best_y,
    plot_fit_width, plot_fit_height,
    GEP_MAX_FITNESS,
    stat_params.total_generations_count + 1,
    2,
    "Best");

  gr_plot_mean = VISUALIZE_AddPlot(
    gr_ctx,
    plot_mean_x,    plot_mean_y,
    plot_fit_width, plot_fit_height,
    GEP_MAX_FITNESS,
    stat_params.total_generations_count + 1,
    2,
    "Mean");

  gr_lb_progress = VISUALIZE_AddLabel(
    gr_ctx,
    label_x, label_y,
    "Progress: 0");
}


/*----------------------------------------------------------------------------*/
//!
static void* gep_run_routine(void *arg)
{
  GEP_Run(gep_ctx, stat_callback);
  graph_to_stop = 1;
  while (!graph_stopped);

  arg = arg;
  return NULL;
}


/*----------------------------------------------------------------------------*/
//!
static gep_short_params_t* createDefaultParams(void)
{
  gep_short_params_t *short_params = malloc(sizeof(*short_params));

  memset(short_params, 0, sizeof(*short_params));

  short_params->population_size          = 50;
  short_params->tree_depth               = 4;
  short_params->genes_count              = 1;
  short_params->generations_count        = 1000;
  short_params->mutations_per_chromosome = 2;
  short_params->ops_preset               = GEP_OPLIST__SHORT;
  short_params->coding_type              = GEP_CODING__BREADTH_FIRST;
  short_params->fitness_type             = GEP_FITNESS__MSE;
  short_params->mse_coefficient          = 1;

  short_params->use_tournament                    = 0;
  short_params->use_replace_worst                 = 0;
  short_params->use_incremental_evolution         = 0;
  short_params->use_selection_probability_density = 0;
  short_params->use_additional_population         = 0;
  short_params->use_dynamic_constants             = 0;

  short_params->input_samples_filename[0]  = '\0';
  short_params->output_samples_filename[0] = '\0';

  return short_params;
}


/*----------------------------------------------------------------------------*/
//!
static int readConfigFile(const char *filename, gep_short_params_t *short_params)
{
  FILE *cfg_file = fopen(filename, "r");

  if (NULL == cfg_file)
  {
    fprintf(stderr, "Can't open config file\n");
    return -1;
  }

  char ops_str[128]          = "";
  char coding_type_str[128]  = "";
  char fitness_type_str[128] = "";

  while (!feof(cfg_file))
  {
    char line[256];
    char token_buf[128];

    if (NULL != fgets(line, sizeof(line), cfg_file))
    {
      CHECK_FIRST_LINE_INT(population_size)
      CHECK_NEXT_LINE_INT(tree_depth)
      CHECK_NEXT_LINE_INT(genes_count)
      CHECK_NEXT_LINE_INT(generations_count)
      CHECK_NEXT_LINE_INT(mutations_per_chromosome)
      CHECK_NEXT_LINE_DBL(mse_coefficient)
      CHECK_NEXT_LINE_INT(use_tournament)
      CHECK_NEXT_LINE_INT(use_replace_worst)
      CHECK_NEXT_LINE_INT(use_incremental_evolution)
      CHECK_NEXT_LINE_INT(use_selection_probability_density)
      CHECK_NEXT_LINE_INT(use_additional_population)
      CHECK_NEXT_LINE_INT(use_dynamic_constants)
      CHECK_NEXT_LINE_STR(ops_preset,   ops_str)
      CHECK_NEXT_LINE_STR(coding_type,  coding_type_str)
      CHECK_NEXT_LINE_STR(fitness_type, fitness_type_str)
    }
  }

  fclose(cfg_file);

  if (0 != strlen(ops_str))
  {
    if (0 == strcmp("short", ops_str))
    {
      short_params->ops_preset = GEP_OPLIST__SHORT;
    }
    else if (0 == strcmp("middle", ops_str))
    {
      short_params->ops_preset = GEP_OPLIST__MIDDLE;
    }
    else if (0 == strcmp("big", ops_str))
    {
      short_params->ops_preset = GEP_OPLIST__BIG;
    }
    else if (0 == strcmp("full", ops_str))
    {
      short_params->ops_preset = GEP_OPLIST__FULL;
    }
    else
    {
      fprintf(stderr, "Unsupported functional preset (%s), set to default\n", ops_str);
    }
  }

  if (0 != strlen(coding_type_str))
  {
    if (0 == strcmp("ferreira", coding_type_str))
    {
      short_params->coding_type = GEP_CODING__BREADTH_FIRST;
    }
    else if (0 == strcmp("prefix", coding_type_str))
    {
      short_params->coding_type = GEP_CODING__PREFIX;
    }
    else if (0 == strcmp("overlapped", coding_type_str))
    {
      short_params->coding_type = GEP_CODING__OVERLAPPED;
    }
    else
    {
      fprintf(stderr, "Unsupported coding type (%s), set to default\n", coding_type_str);
    }
  }

  if (0 != strlen(fitness_type_str))
  {
    if (0 == strcmp("mse", fitness_type_str))
    {
      short_params->fitness_type = GEP_FITNESS__MSE;
    }
    else if (0 == strcmp("partial_mse", fitness_type_str))
    {
      short_params->fitness_type = GEP_FITNESS__PARTIAL_MSE;
    }
    else if (0 == strcmp("r_squared", fitness_type_str))
    {
      short_params->fitness_type = GEP_FITNESS__R_SQUARED;
    }
    else
    {
      fprintf(stderr, "Unsupported fitness type (%s), set to default\n", fitness_type_str);
    }
  }

  return 0;
}


/*----------------------------------------------------------------------------*/
//!
int main(int argc, char *argv[])
{
  if (argc != 5)
  {
    fprintf(stderr, "Usage: gep_test <config_file> <input_samples_file> <output_samples_file> <enable_graphics>\n");
    exit(0);
  }

  gep_short_params_t *short_params = createDefaultParams();

  const char *arg_config_filename         = argv[1];
  const char *arg_input_samples_filename  = argv[2];
  const char *arg_output_samples_filename = argv[3];
  const char *arg_enable_graphics_str     = argv[4];

  readConfigFile(arg_config_filename, short_params);
  strncpy(short_params->input_samples_filename,  arg_input_samples_filename,  sizeof(short_params->input_samples_filename));
  strncpy(short_params->output_samples_filename, arg_output_samples_filename, sizeof(short_params->output_samples_filename));

  sscanf(arg_enable_graphics_str, "%d",  &graph_enabled);

  gep_ctx = GEP_CreateContext(short_params, &stat_params);

  if (signal(SIGUSR1, gep_sigterm_handler) == SIG_ERR)
  {
    fprintf(stderr, "Can't setup SIGUSR1 hanlder\n");
  }

  last_frame_rendered = simulation_started = last_output_printed = getTimeNs();

  if (graph_enabled)
  {
    graph_init(argc, argv);

    pthread_t gep_thread;
    int status = pthread_create(&gep_thread, NULL, gep_run_routine, NULL);
    if (0 != status)
    {
      fprintf(stderr, "Can't create GEP thread!\n");
      exit(-1);
    }

    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
    glutMainLoop();

    pthread_join(gep_thread, NULL);
  }
  else
  {
    GEP_Run(gep_ctx, stat_callback);
  }

  GEP_DestroyContext(gep_ctx);
  free(short_params);

  return 0;
}
