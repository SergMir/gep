#include <visualize.h>

#include <GL/glut.h>
#include <memory.h>
#include <string.h>
#include <pthread.h>

#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 800
#define WORLD_WIDTH 1000
#define WORLD_HEIGHT 1000

#define DARK 1

/*----------------------------------------------------------------------------*/
//!
struct visualize_plot_s
{
  double  x;
  double  y;
  double  width;
  double  height;
  double  max_height;
  int     points_count;
  int     printed_points_count;
  double *points;
  char    title[VISUALIZE_MAX_LABEL_LEN];
  int     dimensions_count;
  int     rotation_alpha;
  visualize_ctx_t
         *ctx; /* Back link */
};

/*----------------------------------------------------------------------------*/
//!
struct visualize_label_s
{
  double  x;
  double  y;
  char    text[VISUALIZE_MAX_LABEL_LEN];
  visualize_ctx_t
         *ctx; /* Back link */
};

/*----------------------------------------------------------------------------*/
//!
struct visualize_ctx_s
{
  int    valid;
  int    gl_window;
  pthread_mutex_t
         mutex;

  int    plots_count;
  visualize_plot_t
         plots[VISUALIZE_MAX_OBJ_CNT];
  double plots_rotation_alpha;

  int    labels_count;
  visualize_label_t
         labels[VISUALIZE_MAX_OBJ_CNT];
};



/*----------------------------------------------------------------------------*/
//!
void visualize_Reshape(int width, int height)
{
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, WORLD_WIDTH, 0, WORLD_HEIGHT, -WORLD_WIDTH, WORLD_WIDTH);
  glMatrixMode(GL_MODELVIEW);
}


/*----------------------------------------------------------------------------*/
//!
visualize_ctx_t* VISUALIZE_CreateContext(void)
{
  visualize_ctx_t *ctx = malloc(sizeof(*ctx));

  glutInitWindowSize(800, 600);
  glutInitWindowPosition(100, 100);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  ctx->gl_window = glutCreateWindow("Gene Expression Programming");

  glutReshapeFunc(visualize_Reshape);
  glClearColor(0, 0, 0, 0);

  glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glEnable(GL_DEPTH_TEST);


  ctx->plots_count  = 0;
  ctx->labels_count = 0;
  memset(ctx->plots,  0, sizeof(ctx->plots));
  memset(ctx->labels, 0, sizeof(ctx->labels));

  pthread_mutex_init(&ctx->mutex, NULL);

  ctx->valid = 1;

  return ctx;
}


/*----------------------------------------------------------------------------*/
//!
void VISUALIZE_DestroyContext(visualize_ctx_t *ctx)
{
  pthread_mutex_lock(&ctx->mutex);

  ctx->valid = 0;

  glutReshapeFunc(NULL);
  glutDestroyWindow(ctx->gl_window);

  int plots_count   = ctx->plots_count;
  ctx->plots_count  = 0;
  ctx->labels_count = 0;

  for (int i = 0; i < plots_count; ++i)
  {
    free(ctx->plots[i].points);
  }
  memset(ctx->plots,  0, sizeof(ctx->plots));
  memset(ctx->labels, 0, sizeof(ctx->labels));

  pthread_mutex_unlock(&ctx->mutex);

  pthread_mutex_destroy(&ctx->mutex);

  free(ctx);
}


/*----------------------------------------------------------------------------*/
//!
void visualize_DrawText(double x, double y, double size, const char *text)
{ 
  glPushMatrix();
  {
  #if defined(DARK)
    glColor3f(1.0f, 1.0f, 1.0f);
  #else
    glColor3f(0, 0, 0);
  #endif

    glTranslatef(x, y, 0);
    glScalef(0.2 * size, 0.2 * size, 0.2 * size);

    for (const char *c = text; *c != '\0'; ++c)
    {
      glutStrokeCharacter(GLUT_STROKE_ROMAN, *c);
    }
  }
  glPopMatrix();
}


/*----------------------------------------------------------------------------*/
//!
void visualize_DrawPlot(visualize_plot_t *plot)
{
  glPushMatrix();

  glTranslated(plot->x, plot->y, 0);

  if (2 == plot->dimensions_count)
  {
    double x_delta = plot->width / plot->points_count;

    glColor3f(0, 1.0f, 0);
    glBegin(GL_LINE_STRIP);
    glVertex2f(0, plot->height);
    glVertex2f(0, 0);
    glVertex2f(plot->width, 0);
    glEnd();

    glColor3f(1.0f, 0, 0);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < plot->printed_points_count; ++i)
    {
      glVertex2f(x_delta * i, (plot->points[i] / plot->max_height) * plot->height);
    }
    glEnd();
  }
  else if (3 == plot->dimensions_count)
  {
    glTranslated(plot->width / 2, plot->height / 2, 0);
    glRotated(-75, 1.0, 0, 0);
    glRotated(plot->rotation_alpha, 0, 0, 1.0);
    glTranslated(- plot->width / 2, - plot->height / 2, 0);

    for (int x = 0; x < (plot->points_count - 1); ++x)
    {
      for (int y = 0; y < (plot->points_count - 1); ++y)
      {
        double *p1 = plot->points + 3 * (y * plot->points_count + x);
        double *p2 = plot->points + 3 * (y * plot->points_count + x + 1);
        double *p3 = plot->points + 3 * ((y + 1) * plot->points_count + x + 1);
        double *p4 = plot->points + 3 * ((y + 1) * plot->points_count + x);

        glColor3f(1.0f, 0, 0);
        glBegin(GL_LINE_LOOP);
        glVertex3f(p1[0], p1[1], (p1[2] / plot->max_height) * plot->height);
        glVertex3f(p2[0], p2[1], (p2[2] / plot->max_height) * plot->height);
        glVertex3f(p3[0], p3[1], (p3[2] / plot->max_height) * plot->height);
        glVertex3f(p4[0], p4[1], (p4[2] / plot->max_height) * plot->height);
        glEnd();

        glColor3f(0, 1.0f, 0);
        glBegin(GL_LINE_LOOP);
        glVertex3f(p1[0], p1[1], 0);
        glVertex3f(p2[0], p2[1], 0);
        glVertex3f(p3[0], p3[1], 0);
        glVertex3f(p4[0], p4[1], 0);
        glEnd();
      }
    }
  }
  glPopMatrix();

  visualize_DrawText(plot->x, plot->y - 20, 0.8, plot->title);
}


/*----------------------------------------------------------------------------*/
//!
void VISUALIZE_Render(visualize_ctx_t *ctx)
{
  pthread_mutex_lock(&ctx->mutex);

  if (ctx->valid)
  {
    glLoadIdentity();
    #if defined(DARK)
      glClearColor(0, 0, 0, 0);
    #else
      glClearColor(1.0f, 1.0f, 1.0f, 0);
    #endif
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    for (int i = 0; i < ctx->plots_count; ++i)
    {
      visualize_plot_t *plot = &ctx->plots[i];
    
      plot->rotation_alpha += 10;
      if (plot->rotation_alpha > 360)
      {
        plot->rotation_alpha = 0;
      }

      visualize_DrawPlot(plot);
    }

    for (int i = 0; i < ctx->labels_count; ++i)
    {
      visualize_label_t *label = &ctx->labels[i];
      visualize_DrawText(label->x, label->y, 0.8, label->text);
    }

    glutSwapBuffers();
  }

  pthread_mutex_unlock(&ctx->mutex);
}


/*----------------------------------------------------------------------------*/
//!
visualize_plot_t* VISUALIZE_AddPlot(visualize_ctx_t *ctx, double x, double y, double width, double height, double max_height, int points_per_dimension, int dimensions_count, const char *title)
{
  visualize_plot_t *plot = NULL;

  pthread_mutex_lock(&ctx->mutex);

  if (ctx->valid)
  {
    plot = &ctx->plots[ctx->plots_count++];

    plot->x                    = x;
    plot->y                    = y;
    plot->width                = width;
    plot->height               = height;
    plot->max_height           = max_height;
    plot->points_count         = points_per_dimension;
    plot->printed_points_count = 0;
    plot->dimensions_count     = dimensions_count;
    plot->rotation_alpha       = 0;
    plot->ctx                  = ctx;
    strncpy(plot->title, title, VISUALIZE_MAX_LABEL_LEN);

    if (2 == plot->dimensions_count)
    {
      plot->points = malloc(sizeof(plot->points[0]) * plot->points_count);
      memset(plot->points, 0, sizeof(plot->points[0]) * plot->points_count);
    }
    else if (3 == plot->dimensions_count)
    {
      plot->points = malloc(3 * sizeof(plot->points[0]) * plot->points_count * plot->points_count);

      for (int y = 0, i = 0; y < points_per_dimension; ++y)
      {
        for (int x = 0; x < points_per_dimension; ++x, ++i)
        {
          double *point = plot->points + 3 * i;
          point[0] = plot->width * (double)x / (double)plot->points_count;
          point[1] = plot->height * (double)y / (double)plot->points_count;
          point[2] = 0;
        }
      }
    }
  }

  pthread_mutex_unlock(&ctx->mutex);

  return plot;
}


/*----------------------------------------------------------------------------*/
//!
visualize_label_t* VISUALIZE_AddLabel(visualize_ctx_t *ctx, double x, double y, const char *text)
{
  visualize_label_t *label = NULL;

  pthread_mutex_lock(&ctx->mutex);

  if (ctx->valid)
  {
    label = &ctx->labels[ctx->labels_count++];

    label->x   = x;
    label->y   = y;
    label->ctx = ctx;
    strncpy(label->text, text, VISUALIZE_MAX_LABEL_LEN);
  }

  pthread_mutex_unlock(&ctx->mutex);

  return label;
}


/*----------------------------------------------------------------------------*/
//!
void VISUALIZE_PlotAdd2dPoint(visualize_plot_t *plot, double y)
{
  visualize_ctx_t *ctx = plot->ctx;

  pthread_mutex_lock(&ctx->mutex);

  if (ctx->valid)
  {
    if (2 == plot->dimensions_count && plot->printed_points_count < plot->points_count)
    {
      plot->points[plot->printed_points_count++] = y;
    }
  }

  pthread_mutex_unlock(&ctx->mutex);
}


/*----------------------------------------------------------------------------*/
//!
void VISUALIZE_PlotSetPoints(visualize_plot_t *plot, const double *points, int points_count)
{
  visualize_ctx_t *ctx = plot->ctx;

  pthread_mutex_lock(&ctx->mutex);

  if (ctx->valid)
  {
    if (2 == plot->dimensions_count)
    {
      memcpy(plot->points, points, sizeof(double) * points_count);
      plot->printed_points_count = points_count;
    }
    else if (3 == plot->dimensions_count)
    {
      for (int i = 0; i < points_count; ++i)
      {
        plot->points[3 * i + 2] = points[i];
      }
    }
  }

  pthread_mutex_unlock(&ctx->mutex);
}


/*----------------------------------------------------------------------------*/
//!
void VISUALIZE_LabelSetText(visualize_label_t *label, const char *text)
{
  visualize_ctx_t *ctx = label->ctx;

  pthread_mutex_lock(&ctx->mutex);

  if (ctx->valid)
  {
    strncpy(label->text, text, VISUALIZE_MAX_LABEL_LEN);
  }

  pthread_mutex_unlock(&ctx->mutex);
}
