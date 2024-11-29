#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(gettree)(int *iTree, int *nr, int *nc, double *nodes, double *Func, double *mean, double *Prob);
extern void F77_NAME(predictsurvtree)(int *n, int *np, double *xt, int *nCat, int *nt, int *nNodes, double *tFunc, double *mean, double *Prob, int *nCols, double *tnodes, double *predFunc, double *predMean, double *predProb);
extern void F77_NAME(setupbasics)(double *t_surv_tp, double *t_end_tp, int *t_nt, int *t_nt_death, double *t_dt, double *t_dt_death, double *t_rs, int *t_ERT, int *t_uniformSplit, int *t_nodeSizeEnd, int *t_nodeSizeSurv, int *t_minEventEnd, int *t_minEventSurv, int *t_rule, int *t_sIndex, double *t_sFraction, double *t_stratifiedSplit, int *t_replace);
extern void F77_NAME(setupinners)(int *t_n, int *t_n_surv, int *t_idvec, int *t_person_ind, int *t_np, double *t_x, double *t_pr, double *t_pr2, double *t_pr2surv, double *t_prsurv, double *t_ord_causeind, double *t_ord_response, int *t_delta, int *t_delta_m, int *t_mTry, int *t_nCat, int *t_sampleSize, int *t_nTree, int *t_nrNodes);
extern void F77_NAME(survtree)(double *tSurvFunc, double *mean, double *survProb);
extern void F77_NAME(endpointtree)(double *tFunc, double *mean, double *Prob);
extern void F77_NAME(treedim)(int *iTree, int *nr, int *nc);
extern void F77_NAME(temporary)(double *s, double *vs, double *z);
extern void F77_NAME(crstm)(double *y, int *m, int *ig, int *ist, int *no, double *rho, int *nst, int *ng, double *s, double *vs, double *ys, int *ms, int *igs, double *v, double *st, double *vt, double *wk, int *iwk, double *z);
extern void F77_NAME(example_sort)(int *n);

static const R_FortranMethodDef FortranEntries[] = {
    {"gettree",         (DL_FUNC) &F77_NAME(gettree),          7},
    {"predictsurvtree", (DL_FUNC) &F77_NAME(predictsurvtree), 14},
    {"setupbasics",     (DL_FUNC) &F77_NAME(setupbasics),     18},
    {"setupinners",     (DL_FUNC) &F77_NAME(setupinners),     19},
    {"survtree",        (DL_FUNC) &F77_NAME(survtree),         3},
    {"endpointtree",    (DL_FUNC) &F77_NAME(endpointtree),         3},
    {"treedim",         (DL_FUNC) &F77_NAME(treedim),          3},
    {"temporary",       (DL_FUNC) &F77_NAME(temporary),        3},
    {"crstm",           (DL_FUNC) &F77_NAME(crstm),            18},
    {"example_sort",    (DL_FUNC) &F77_NAME(example_sort),     1},
    {NULL, NULL, 0}
};

void R_init_itrSurv(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

