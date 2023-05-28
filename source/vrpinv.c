/*  S-Function for Simulink Real Matrix Pseudoinversion ***************************************/
/*  Giampiero Campa 27-August-00 **************************************************************/

#define S_FUNCTION_NAME vrpinv
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(a,b)  ((a) > (b)  ? (a) : (b))
#define MIN(a,b)  ((a) < (b)  ? (a) : (b))
#define INV(a,b)  ((a) > (b)  ? (1/a) : (0.0))

real_T pythag2(real_T a, real_T b)
{
	real_T at=fabs(a), bt=fabs(b), ct;
	return ((at=fabs(a)) > (bt=fabs(b)) ? \
	(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0));
}

void svdcmp2(real_T* rwv,int_T m,int_T n)
{
	int_T flag,i,its,j,jj,k,l,nm;
	real_T c,f,h,s,x,y,z;
	real_T anorm=0.0,g=0.0,scale=0.0;

	for (i=1;i<=n;i++) {
		l=i+1;
		rwv[n*(m+n+1)-1+i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(rwv[-1+k-m+m*i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					rwv[-1+k-m+m*i] /= scale;
					s += rwv[-1+k-m+m*i]*rwv[-1+k-m+m*i];
				}
				f=rwv[-1+i-m+m*i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				rwv[-1+i-m+m*i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += rwv[-1+k-m+m*i]*rwv[-1+k-m+m*j];
						f=s/h;
						for (k=i;k<=m;k++) rwv[-1+k-m+m*j] += f*rwv[-1+k-m+m*i];
					}
				}
				for (k=i;k<=m;k++) rwv[-1+k-m+m*i] *= scale;
			}
		}
		rwv[m*n-1+i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(rwv[-1+i-m+m*k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					rwv[-1+i-m+m*k] /= scale;
					s += rwv[-1+i-m+m*k]*rwv[-1+i-m+m*k];
				}
				f=rwv[-1+i-m+m*l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				rwv[-1+i-m+m*l]=f-g;
				for (k=l;k<=n;k++) rwv[n*(m+n+1)-1+k]=rwv[-1+i-m+m*k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += rwv[-1+j-m+m*k]*rwv[-1+i-m+m*k];
						for (k=l;k<=n;k++) rwv[-1+j-m+m*k] += s*rwv[n*(m+n+1)-1+k];
					}
				}
				for (k=l;k<=n;k++) rwv[-1+i-m+m*k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(rwv[m*n-1+i])+fabs(rwv[n*(m+n+1)-1+i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					rwv[n*(m+1)-1+j-n+n*i]=(rwv[-1+i-m+m*j]/rwv[-1+i-m+m*l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += rwv[-1+i-m+m*k]*rwv[n*(m+1)-1+k-n+n*j];
					for (k=l;k<=n;k++) rwv[n*(m+1)-1+k-n+n*j] += s*rwv[n*(m+1)-1+k-n+n*i];
				}
			}
			for (j=l;j<=n;j++) rwv[n*(m+1)-1+i-n+n*j]=rwv[n*(m+1)-1+j-n+n*i]=0.0;
		}
		rwv[n*(m+1)-1+i-n+n*i]=1.0;
		g=rwv[n*(m+n+1)-1+i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=rwv[m*n-1+i];
		if (i < n)
			for (j=l;j<=n;j++) rwv[-1+i-m+m*j]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += rwv[-1+k-m+m*i]*rwv[-1+k-m+m*j];
					f=(s/rwv[-1+i-m+m*i])*g;
					for (k=i;k<=m;k++) rwv[-1+k-m+m*j] += f*rwv[-1+k-m+m*i];
				}
			}
			for (j=i;j<=m;j++) rwv[-1+j-m+m*i] *= g;
		} else {
			for (j=i;j<=m;j++) rwv[-1+j-m+m*i]=0.0;
		}
		++rwv[-1+i-m+m*i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((real_T)(fabs(rwv[n*(m+n+1)-1+l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((real_T)(fabs(rwv[m*n-1+nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rwv[n*(m+n+1)-1+i];
					rwv[n*(m+n+1)-1+i]=c*rwv[n*(m+n+1)-1+i];
					if ((real_T)(fabs(f)+anorm) == anorm) break;
					g=rwv[m*n-1+i];
					h=pythag2(f,g);
					rwv[m*n-1+i]=h;
					h=1.0/h;
					c=g*h;
					s=(-f*h);
					for (j=1;j<=m;j++) {
						y=rwv[-1+j-m+m*nm];
						z=rwv[-1+j-m+m*i];
						rwv[-1+j-m+m*nm]=y*c+z*s;
						rwv[-1+j-m+m*i]=z*c-y*s;
					}
				}
			}
			z=rwv[m*n-1+k];
			if (l == k) {
				if (z < 0.0) {
					rwv[m*n-1+k] = -z;
					for (j=1;j<=n;j++) rwv[n*(m+1)-1+j-n+n*k]=(-rwv[n*(m+1)-1+j-n+n*k]);
				}
				break;
			}
			if (its == 50) 
			{
				printf("No convergence in 50 SVDCMP2 iterations \n");
				return;
			}

			x=rwv[m*n-1+l];
			nm=k-1;
			y=rwv[m*n-1+nm];
			g=rwv[n*(m+n+1)-1+nm];
			h=rwv[n*(m+n+1)-1+k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag2(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rwv[n*(m+n+1)-1+i];
				y=rwv[m*n-1+i];
				h=s*g;
				g=c*g;
				z=pythag2(f,h);
				rwv[n*(m+n+1)-1+j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=rwv[n*(m+1)-1+jj-n+n*j];
					z=rwv[n*(m+1)-1+jj-n+n*i];
					rwv[n*(m+1)-1+jj-n+n*j]=x*c+z*s;
					rwv[n*(m+1)-1+jj-n+n*i]=z*c-x*s;
				}
				z=pythag2(f,h);
				rwv[m*n-1+j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=rwv[-1+jj-m+m*j];
					z=rwv[-1+jj-m+m*i];
					rwv[-1+jj-m+m*j]=y*c+z*s;
					rwv[-1+jj-m+m*i]=z*c-y*s;
				}
			}
			rwv[n*(m+n+1)-1+l]=0.0;
			rwv[n*(m+n+1)-1+k]=f;
			rwv[m*n-1+k]=x;
		}
	}
}

/* mdlCheckParameters, check parameters, this routine is called later from mdlInitializeSizes */
#define MDL_CHECK_PARAMETERS
static void mdlCheckParameters(SimStruct *S)
{
    /* Basic check : All parameters must be integer positive vectors                             */
    real_T *pr;                            

    int_T  i, el, nEls;
    for (i = 0; i < 1; i++) {
        if (mxIsEmpty(    ssGetSFcnParam(S,i)) || mxIsSparse(   ssGetSFcnParam(S,i)) ||
            mxIsComplex(  ssGetSFcnParam(S,i)) || !mxIsNumeric( ssGetSFcnParam(S,i))  )
                  { ssSetErrorStatus(S,"Parameters must be real finite vectors"); return; } 
        pr   = mxGetPr(ssGetSFcnParam(S,i));
        nEls = mxGetNumberOfElements(ssGetSFcnParam(S,i));
        for (el = 0; el < nEls; el++) {
            if (!mxIsFinite(pr[el])) 
                  { ssSetErrorStatus(S,"Parameters must be real finite vectors"); return; }
        }
    }

    /* Check number of elements in parameter: [no ni]                                         */
    if ( mxGetNumberOfElements(ssGetSFcnParam(S,0)) != 2 )
    { ssSetErrorStatus(S,"The parameter must be a 2 elements vector"); return; }

    /* get the basic parameters and check them                                                */
    pr=mxGetPr(ssGetSFcnParam(S,0));
    if ( (pr[0] < 1) | (pr[1] < 1) )
    { ssSetErrorStatus(S,"Dimensions must be greater than zero"); return; }

    /* Check number of elements in parameter: tol                                             */
    if ( mxGetNumberOfElements(ssGetSFcnParam(S,1)) != 1 )
    { ssSetErrorStatus(S,"The parameter must be a real number"); return; }

    /* get the basic parameters and check them                                                */
    pr=mxGetPr(ssGetSFcnParam(S,1));
    if ( pr[0] < 0 )
    { ssSetErrorStatus(S,"Tolerance must be greater than zero"); return; }

}

/* mdlInitializeSizes - initialize the sizes array ********************************************/
static void mdlInitializeSizes(SimStruct *S)
{
    real_T *n;                            
    int_T erw;

    ssSetNumSFcnParams(S,2);                          /* number of expected parameters        */

    /* Check the number of parameters and then calls mdlCheckParameters to see if they are ok */
    if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S))
    { mdlCheckParameters(S); if (ssGetErrorStatus(S) != NULL) return; } else return;
    n=mxGetPr(ssGetSFcnParam(S,0));

    ssSetNumContStates(S,0);                          /* number of continuous states          */
    ssSetNumDiscStates(S,0);                          /* number of discrete states            */

    if (!ssSetNumInputPorts(S,1)) return;             /* number of input ports                */
    ssSetInputPortWidth(S,0,(int_T)(n[0]*n[1]));      /* first input port width               */
    ssSetInputPortDirectFeedThrough(S,0,1);           /* first port direct feedthrough flag   */

    if (!ssSetNumOutputPorts(S,1)) return;            /* number of output ports               */
    ssSetOutputPortWidth(S,0,(int_T)(n[1]*n[0]));     /* first output port width              */
   
    ssSetNumSampleTimes(S,0);                         /* number of sample times               */

    /* number real_T work vector elements													  */
	erw=(int_T)((n[0]+n[1]+2)*MIN(n[1],n[0]));
    ssSetNumRWork(S,erw+1);

    ssSetNumIWork(S,2);                               /* number int_T work vector elements    */
    ssSetNumPWork(S,0);                               /* number ptr work vector elements      */
    ssSetNumModes(S,0);                               /* number mode work vector elements     */
    ssSetNumNonsampledZCs(S,0);                       /* number of nonsampled zero crossing   */
}

/* mdlInitializeSampleTimes - initialize the sample times array *******************************/
static void mdlInitializeSampleTimes(SimStruct *S)
{
    /* Set things up to run with inherited sample time                                        */
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0);
}

/* mdlStart - initialize work vectors *********************************************************/
#define MDL_START
static void mdlStart(SimStruct *S)
{
real_T           *prv = mxGetPr(ssGetSFcnParam(S,0));
real_T   *t, *n, *rwv = ssGetRWork(S);
int_T    i, erw, *iwv = ssGetIWork(S);

for (i=0;i<2;i++) iwv[i]=(int_T)(prv[i]);

n=mxGetPr(ssGetSFcnParam(S,0));
erw=(int_T)((n[0]+n[1]+2)*MIN(n[1],n[0]));

for (i=0;i<erw;i++) rwv[i]=(0.0);

t=mxGetPr(ssGetSFcnParam(S,1));
rwv[erw]=t[0];
}

/* mdlOutputs - compute the outputs ***********************************************************/
static void mdlOutputs(SimStruct *S, int_T tid)
{
int_T         i, j, k, *n  = ssGetIWork(S);
real_T                 *y  = ssGetOutputPortRealSignal(S,0);

real_T			   t, *rwv = ssGetRWork(S);

InputRealPtrsType      up  = ssGetInputPortRealSignalPtrs(S,0);

if (n[0]>=n[1])
{	
	// A=reshape(1:6,3,2);[U,S,V]=svd(A,0); A-U*S*V',pinv(A)-V*inv(S)*U'

	for(j = 0; j < n[1]*n[0]; j++) rwv[j] = (*up[j]);

	svdcmp2(rwv,n[0],n[1]);

	t=n[0]*rwv[n[1]*n[0]]*2.220446049250313e-016*rwv[(n[0]+n[1]+2)*n[1]];

	for(i = 0; i < n[1]; i++)
		for(j = 0; j < n[0]; j++)
			for( y[i+j*n[1]] = 0, k = 0; k < n[1]; k++)
				 y[i+j*n[1]] += rwv[n[1]*(n[0]+1)+i+k*n[1]]*INV(rwv[n[1]*n[0]+k],t)*rwv[j+k*n[0]];
}
else
{
	// A=reshape(1:6,2,3);[V,S,U]=svd(A',0); A-U*S*V',pinv(A)-V*inv(S)*U'

	for(i = 0; i < n[1]; i++)
		for(j = 0; j < n[0]; j++)
			rwv[i+j*n[1]] = (*up[j+i*n[0]]);

	svdcmp2(rwv,n[1],n[0]);

	t=n[1]*rwv[n[1]*n[0]]*2.220446049250313e-016*rwv[(n[0]+n[1]+2)*n[0]];

	for(i = 0; i < n[1]; i++)
		for(j = 0; j < n[0]; j++)
			for( y[i+j*n[1]] = 0, k = 0; k < n[0]; k++)
				 y[i+j*n[1]] += rwv[i+k*n[1]]*INV(rwv[n[1]*n[0]+k],t)*rwv[n[0]*(n[1]+1)+j+k*n[0]];
}

}

/* mdlTerminate - called when the simulation is terminated ***********************************/
static void mdlTerminate(SimStruct *S) {}

/* Trailer information to set everything up for simulink usage *******************************/
#ifdef  MATLAB_MEX_FILE                      /* Is this file being compiled as a MEX-file?   */
#include "simulink.c"                        /* MEX-file interface mechanism                 */
#else
#include "cg_sfun.h"                         /* Code generation registration function        */
#endif

#undef SIGN
#undef MAX
#undef MIN
#undef INV

