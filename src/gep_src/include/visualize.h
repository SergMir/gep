#ifndef _H_VISUALIZE_
#define _H_VISUALIZE_

#define VISUALIZE_MAX_OBJ_CNT   10
#define VISUALIZE_MAX_LABEL_LEN 256

typedef struct visualize_ctx_s   visualize_ctx_t;
typedef struct visualize_plot_s  visualize_plot_t;
typedef struct visualize_label_s visualize_label_t;

visualize_ctx_t* VISUALIZE_CreateContext(void);
void VISUALIZE_DestroyContext(visualize_ctx_t *ctx_handle);
void VISUALIZE_Render        (visualize_ctx_t *ctx_handle);

visualize_plot_t*  VISUALIZE_AddPlot (visualize_ctx_t *ctx_handle, double x, double y, double width, double height, double max_height, int points_per_dimension, int dimensions_count, const char *title);
visualize_label_t* VISUALIZE_AddLabel(visualize_ctx_t *ctx_handle, double x, double y, const char *text);

void VISUALIZE_PlotAdd2dPoint(visualize_plot_t *plot, double y);
void VISUALIZE_PlotSetPoints (visualize_plot_t *plot, const double *points, int points_count);

void VISUALIZE_LabelSetText(visualize_label_t *label, const char *text);

#endif /* _H_VISUALIZE_ */
