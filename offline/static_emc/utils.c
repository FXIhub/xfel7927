
// this is called once for each class and frame sub-set
// one worker per frame

// LR[d]   = sum_i (K[i, d] log(T[i, d]) - T[i, d])
// T[i, d] = w[d] W[i] + sum_l b[d, l] B[l, i]

kernel void calculate_LR_T_dt (
global float *LR,  
global unsigned char *K,  
global float *w, 
global float *W,
global float *b, 
global float *B, 
const float beta, 
const int L, 
const int I, 
const int D)
{   
int d = get_global_id(0);

float T;
float LR_local = 0. ;

float wd = w[d];

int i, l;

for (i=0; i<I; i++){
    // calculate T
    T = wd * W[i] ;
    
    // add background
    for (l=0; l<L; l++) 
        T += b[d * L + l] * B[I*l + i];
    
    LR_local += K[D*i + d] * log(T) - T;
}

LR_local *= beta ;

LR[d] = LR_local ;
}


float take_step(
float x, 
float g,
float c,
bool *mask,
float min_val,
float tol, 
bool *done)
{

float slope = g * g ;
float step, step_max ;

// if the gradient is zero and x is not at the boundary 
// then we are done
if ((slope < tol*tol) || ((g < 0.) && (mask == 0))) {
    //if v: print('gradient within tolerance, done!')
    *done = 1;
    return x;
}

// or the gradient is negative and x is at the boundary
// then we are done
if ((g < 0.) && (mask == 0)) {
    //if v: print('gradient negative at boundary, done')
    *done = 1;
    return x;
}

if ((mask == 0) && (g > 0.)) {   
    //if v: print('positive gradient at boundary, unblocking')
    *mask = 1 ;
}

// check for non-negative curvature
if (c < 0.) 
    step = - slope / c ;
else 
    //if v : print('Warning, non-negative curvature, taking small step in grad')
    step = 1e-5 / sqrt(slope) ;

// get maximum step size before hitting min_val
// negative gradient of unmasked x values
if (g < 0.) 
    step_max = (min_val - x) / g ;
else 
    step_max = MAXFLOAT ;


// take step 
x += fmin(step, step_max) * g ;

// block x if it hits boundary
if (step_max <= step) {
    //if v: print(f'hit lower bound, blocking')
    *mask = 0 ;
    x = min_val ;
}
    
return x ;

}


// g    = gw0 + sum_j  PWK[j] / T[j]
// c    = g^2 sum_j PW2K[j] / T[j]^2
// T[j] = w W[j] + back[j]
kernel void update_w_old (
global float *W,  
global float *PWK,  
global float *PW2K,  
global float *back,  
const float gw0,  
global float *wd,  
const float min_val,
const int max_iters,
const float tol,
const int J,
const int d)
{   
float w = wd[d];

float g, c, T, wold, dw;

int iters, j, k;

bool m, done;

if (w <= min_val) 
   m = 0;
else 
   m = 1;

done = 0;

// improve precision
float temp_sum_g[1024];
float temp_sum_c[1024];
int step = (int)ceil((float)J / (float)1024) ;

for (iters = 0; iters < max_iters; iters++) {
    wold = w;
    g = gw0;
    c = 0.;
    // calculate gradient and curvature
    for (k = 0; k < 1024; k++) {
        temp_sum_g[k] = 0.;
        temp_sum_c[k] = 0.;
    }   
    
    for (k = 0; k < 1024; k++) {
    for (j = step * k; j < min(step * (k+1), J); j++) {
        T  = w * W[j] + back[j];
        temp_sum_g[k] += PWK[j] / T ;
        temp_sum_c[k] -= PW2K[j] / (T * T) ;
    }}
    
    for (k = 0; k < 1024; k++) {
        g += temp_sum_g[k];
        c += temp_sum_c[k];
    }
    
    c *= g * g;
    
    // take step
    w = take_step(w, g, c, &m, min_val, tol, &done);
    
    dw = fabs(w - wold);
    
    //printf("w g0 g c w wold w-wold %f %f %f %f %f %f\n", w, gw0, g, c, wold, dw);
    
    // 1e-6 is the limit for floating point
    if (dw < 1e-4) 
        break ;
    
    if (done)
        break ;
}

//printf("%f %f %d %d %d\n", w, g, d, J, iters);

wd[d] = w;
}


kernel void calculate_gW(
global float *P,  
global float *w,  
global int *ds,  
global int *Ds,  
global float *g_t,  
const int C,
const int D
)
{

int t = get_global_id(0);

int d, n;

// calculate gt 
float c = 0. ;
for (n = 0; n < Ds[t]; n++){ 
    d  = ds[D * t + n] ;
    c += P[C * d + t] * w[d];
}

g_t[t] = c;
}


// background[d, i] = sum_l b[l, d] B0[l, i] / w[d]
kernel void calculate_background (
global float *B,  
global float *b,  
global float *w,  
global float *background,  
const int I,
const int L,
const int D)
{   

int d = get_global_id(0);
int i = get_global_id(1);

float t = 0.;
float wd = w[d];

for (int l=0; l<L; l++) 
    t += b[l * D + d] * B[l * I + i] ;

background[d * I + i] = t / wd;
}



// g[i] =   sum_d P[d] K[d, i] / (W[i] + sum_l b[d, l] B[l, i] / w[d]) - g0
// c[i] = - sum_d P[d] K[d, i] / (W[i] + sum_l b[d, l] B[l, i] / w[d])**2 - g0

kernel void update_W(
global float *P,  
global unsigned char *K,  
global float *b,  
global float *B,  
global float *w,  
global float *W,  
const float c,
const float minval,
const int I,
const int L,
const int D)
{   

int i = get_global_id(0);

int d, iters;
float f, g, T, PK, u, v, xp;

float x  = W[i];

// calculate xmax = sum_d P[d] K[d, i] / sum_d w[d] P[d]
float maxval = 0.;
for (d = 0; d < D; d++){ 
    maxval += P[d] * K[I*d + i] ;
}
maxval /= c ;

x = clamp(x, minval, maxval);

// optimisation loop 3x 
for (iters = 0; iters < 5; iters++){
    // calculate f and g in this notation
    f = 0.;
    g = 0.;
    for (int d = 0; d < D; d++){ 
        T  = 0. ;
        for (int l = 0; l < L; l++)
            T += b[L*d + l] * B[I*l + i] / w[d] ;
        T += x; 
        
        T  = max(T, minval) ;
        PK = P[d] * K[I*d + i] ;
        f += PK / T ;
        g -= PK / (T * T) ;
    }
    
    u  = - f * f / g ;
    v  = - f / g - x ;
    xp = u / c - v ;
    
    x = clamp(xp, minval, maxval) ;
}

W[i] = x;
}

// g[t, i] =   sum_d P[d, t] K[d, i] / (W[t, i] + B[d, i] / w[d]) - g0[t]
// c[t, i] = - sum_d P[d, t] K[d, i] / (W[t, i] + B[d, i] / w[d])**2 - g0[t]
kernel void update_W_old(
global float *P,  
global unsigned char *K,  
global int *ds,  
global int *Ds,  
global float *g_t,  
global float *B,  
global float *W,  
global float *w,  
const float minval,
const int I,
const int C,
const int D)
{   

int t = get_global_id(0);
int i = get_global_id(1);

int d, n, iters;
float f, g, T, PK, u, v, xp;

float x      = W[I * t + i];

float c = g_t[t] ;

// calculate xmax
float maxval = 0.;
for (n = 0; n < Ds[t]; n++){ 
    d  = ds[D * t + n] ;
    maxval += P[C * d + t] * K[I * d + i] / c ;
}

x = clamp(x, minval, maxval);

// optimisation loop 3x 
for (iters = 0; iters < 3; iters++){
    // calculate f and g in this notation
    f = 0.;
    g = 0.;
    for (int n = 0; n < Ds[t]; n++){ 
        // get the n'th relevant frame index for this class
        d = ds[D * t + n] ;

        T  = x + B[I * d + i] ;
        T  = max(T, minval);
        PK = P[C * d + t] * K[I * d + i] ;
        f += PK / T ;
        g -= PK / (T * T) ;
    }
    
    u  = - f * f / g ;
    v  = - f / g - x ;
    xp = u / c - v ;
    
    x = clamp(xp, minval, maxval) ;
}

W[I * t + i] = x;

}


// g0[l] = sum_t (sum_d P[d, t] b[d, l])
//         sum_d b[d, l] (sum_t P[d, t])
//         sum_d b[d, l] 

kernel void calculate_gB(
global float *b,  
global float *g_l,  
const int L,
const int D
)
{

int l = get_global_id(0);

int d, n;

// calculate gt 
float c = 0. ;
for (d = 0; d < D; d++){ 
    c += b[L*d + l];
}

g_l[l] = c;
}

// grad[i] = sum_dt P[d, t] K[d, i] / ( B[i] + (w[d] W[t, i] + sum_l'neql b[d, l'] B[l', i]) / b[d]) - g0[l]
// xmax[i] = sum_dt P[d, t] K[d, i] / g0[l]
//         = sum_d K[d, i] / g0[l]

kernel void update_B(
global float *P,  
global unsigned char *K,  
global int *ds,  
global int *Ds,  
global float *g_l,  
global float *B,  
global float *b,  
global float *W,  
global float *w,  
const float minval,
const int L,
const int I,
const int C,
const int D
){

int l = get_global_id(0);
int i = get_global_id(1);

int d, n, t, lp, iters;
float f, g, T, PK, u, v, xp;

float x      = B[I * l + i];

float c = g_l[t] ;

// calculate xmax 
float maxval = 0.;
for (d = 0; d < D; d++){ 
    maxval += K[I * d + i] / c ;
}

x = clamp(x, minval, maxval);

// optimisation loop 3x 
for (iters = 0; iters < 3; iters++){
    // calculate f and g in this notation
    f = 0.;
    g = 0.;
    // grad[l, i] =  sum_t sum_d K[d, i] P[d, t] / ( B[i] + (w[d] W[t, i] + sum_l'neql b[d, l'] B[l', i]) / b[d]) - g0[l]
    for (t = 0; t < C; t++){ 
    for (n = 0; n < Ds[t]; n++){ 
        // get the n'th relevant frame index for this class
        d = ds[D * t + n] ;
        T = 0.;
        for (lp = 0; lp < L; lp++) {
            if (lp != l) 
                T += b[L*d + lp] * B[I * lp + i] ;
        }
         
        T  += w[d] * W[I * t + i] ;
        T = x + T / b[L * d + l]  ;
        T  = max(T, minval);
        PK = P[C * d + t] * K[I * d + i];
        f += PK / T ;
        g -= PK / (T * T) ;
    }}
    //printf("%.2e %.2e %.2e\n", x, c, f);
    
    u  = - f * f / g ;
    v  = - f / g - x ;
    xp = u / c - v ;
    
    x = clamp(xp, minval, maxval) ;
}

B[I * l + i] = x;

}


// grad[l, i] = sum_dt P[d, t] K[d, i] / ( B[i] + (w[d] W[t, i] + sum_l'neql b[d, l'] B[l', i]) / b[d]) - g0[l]
// xmax[l, i] = sum_dt P[d, t] K[d, i] / g0[l]
//            = sum_d K[d, i] / g0[l]

kernel void update_B_old(
global float *P,  
global unsigned char *K,  
global int *ds,  
global int *Ds,  
global float *g_l,  
global float *B,  
global float *b,  
global float *W,  
global float *w,  
const float minval,
const int L,
const int I,
const int C,
const int D
){

int l = get_global_id(0);
int i = get_global_id(1);

int d, n, t, lp, iters;
float f, g, T, PK, u, v, xp;

float x      = B[I * l + i];

float c = g_l[t] ;

// calculate xmax 
float maxval = 0.;
for (d = 0; d < D; d++){ 
    maxval += K[I * d + i] / c ;
}

x = clamp(x, minval, maxval);

// optimisation loop 3x 
for (iters = 0; iters < 3; iters++){
    // calculate f and g in this notation
    f = 0.;
    g = 0.;
    // grad[l, i] =  sum_t sum_d K[d, i] P[d, t] / ( B[i] + (w[d] W[t, i] + sum_l'neql b[d, l'] B[l', i]) / b[d]) - g0[l]
    for (t = 0; t < C; t++){ 
    for (n = 0; n < Ds[t]; n++){ 
        // get the n'th relevant frame index for this class
        d = ds[D * t + n] ;
        T = 0.;
        for (lp = 0; lp < L; lp++) {
            if (lp != l) 
                T += b[L*d + lp] * B[I * lp + i] ;
        }
         
        T  += w[d] * W[I * t + i] ;
        T = x + T / b[L * d + l]  ;
        T  = max(T, minval);
        PK = P[C * d + t] * K[I * d + i];
        f += PK / T ;
        g -= PK / (T * T) ;
    }}
    //printf("%.2e %.2e %.2e\n", x, c, f);
    
    u  = - f * f / g ;
    v  = - f / g - x ;
    xp = u / c - v ;
    
    x = clamp(xp, minval, maxval) ;
}

B[I * l + i] = x;

}


// Wt[t] = sum_i W[t, i]
kernel void calculate_Wt(
global float *W,
global float *Wt,
const int I
){

int t = get_global_id(0);

float c = 0.;
for (int i = 0; i < I; i++)
    c += W[I*t + i] ;

Wt[t] = c ;
}

// Wt[t] = sum_i W[i, t]
kernel void calculate_Wt2(
global float *W,
global float *Wt,
const int I,
const int C
){

int t = get_global_id(0);

float c = 0.;
for (int i = 0; i < I; i++)
    c += W[C*i + t] ;

Wt[t] = c ;
}


// gw[d] = sum_t P[d, t] Wt[t]
kernel void calculate_gw(
global float *P,
global float *Wt,
global float *gw,
const int C
){

int d = get_global_id(0);

float c = 0.;
for (int t = 0; t < C; t++)
    c += P[C*d + t] * Wt[t] ;

gw[d] = c ;
}

// background[d, i] = sum_l b[d, l] B[l, i]
kernel void calculate_background_di(
global float *B,
global float *b,
global float *background,
const int L,
const int I
){

int d = get_global_id(0);
int i = get_global_id(1);

float c = 0.;
for (int l = 0; l < L; l++)
    c += b[L*d + l] * B[I*l + i] ;

background[I*d + i] = c ;
}

// background[i, d] = sum_l b[d, l] B[i, l]
kernel void calculate_background_id(
global float *B,
global float *b,
global float *background,
const int L,
const int I,
const int D
){

int d = get_global_id(0);
int i = get_global_id(1);

float c = 0.;
for (int l = 0; l < L; l++)
    //c += b[L*d + l] * B[I*l + i] ;
    c += b[L*d + l] * B[L*i + l] ;

background[D*i + d] = c ;
}

// sub-optimal since we are summing over fast scan
//------------------------------------------------
// grad[d] = sum_t P[d, t] sum_i K[d, i] / (w[d] + sum_l b[d, l] B[l, i] / W[t, i]) - g0[d]
// xmax[d] = (sum_t P[d, t]) (sum_i K[d, i]) / g0[d]
//         = (sum_i K[d, i]) / g0[d]


kernel void update_w(
global float *P,
global unsigned char *K,
global float *background,
global float *W,
global float *w,
global int *Ts,
global int *ts,
global float *gw,
const float minval,
const int I,
const int C,
const int D
){
int d = get_global_id(0);

int n, i, t, iters;
float xp, u, v, f, g, PK, T;

float c = gw[d] ;
float x = w[d];

// calculate xmax[d]
// could precalculate sum_i K[d, i]
float maxval = 0.;
for (i = 0; i < I; i++)
    maxval += K[D*i + d];
maxval /= c ;


x = clamp(x, minval, maxval);

// optimisation loop 3x 
for (iters = 0; iters < 5; iters++){
    // calculate f and g in this notation
    f = 0.;
    g = 0.;
    // ts[d, n] = n'th significant class id (t)
    // grad[d] = sum_t P[d, t] sum_i K[i, d] / (w[d] + sum_l b[d, l] B[i, l] / W[i, t]) - g0[d]
    // T       = w[d] + sum_l b[d, l] B[i, l] / W[i, t]
    for (n = 0; n < Ts[d]; n++){ 
        t = ts[C * d + n] ;
        for (i = 0; i < I; i++){ 
            //T  = x + background[I*d + i] / W[I*t + i];
            T  = x + background[D*i + d] / W[C*i + t];
            T  = max(T, minval);
            //PK = P[C*d + t] * K[I*d + i];
            PK = P[C*d + t] * K[D*i + d];
            f += PK / T ;
            g -= PK / (T * T) ;
    }}
    //printf("%.2e %.2e %.2e\n", x, c, f);
    
    u  = - f * f / g ;
    v  = - f / g - x ;
    xp = u / c - v ;
    
    x = clamp(xp, minval, maxval) ;
}

w[d] = x;
}



// gb[l] = sum_i B[i, l]
kernel void calculate_gb(
global float *B,
global float *gb,
const int I,
const int L
){

int l = get_global_id(0);

float c=0.;
for (int i=0; i < I; i++)
    //c += B[I*l + i];
    c += B[L*i + l];

gb[l] = c;
}

// gb[d, l] = sum_t P[d, t] sum_i K[d, i] / (b[d, l] + (w[d] W[t, i] + sum_l'!=l b[d, l'] B[l', i])/B[l, i]) - sum_i B[l, i]
// xmax[d, l] = (sum_t P[d, t]) (sum_i K[d, i]) / gw[l]
//            = sum_i K[d, i] / gw[l]
kernel void update_b(
global float *P,
global unsigned char *K,
global float *background,
global float *W,
global float *w,
global float *B,
global float *b,
global int *Ts,
global int *ts,
global float *gb,
const float minval,
const int I,
const int C,
const int D,
const int L
){
int d = get_global_id(0);
int l = get_global_id(1);

int n, i, t, iters, lp;
float xp, u, v, f, g, PK, T;

float c = gb[l] ;
float x = b[L*d + l];

// calculate xmax[d]
// should precalculate sum_i K[i, d]
float maxval = 0.;
for (i = 0; i < I; i++)
    maxval += K[D*i + d];
maxval /= c ;


x = clamp(x, minval, maxval);

// optimisation loop 3x 
for (iters = 0; iters < 3; iters++){
    // calculate f and g in this notation
    f = 0.;
    g = 0.;
    // ts[d, n] = n'th significant class id (t)
    // gb[d, l] = sum_t P[d, t] sum_i K[d, i] / (b[d, l] + (w[d] W[t, i] + sum_l'!=l b[d, l'] B[l', i])/B[l, i]) - sum_i B[l, i]
    for (n = 0; n < Ts[d]; n++){ 
        t = ts[C * d + n] ;
        for (i = 0; i < I; i++){ 
            T = 0.;
            for (lp = 0; lp < L; lp++) {
                if (lp != l) 
                    //T += b[L*d + lp] * B[I * lp + i] ;
                    T += b[L*d + lp] * B[L * i + lp] ;
            }
             
            //T  += w[d] * W[I * t + i] ;
            T  += w[d] * W[C * i + t] ;
            //T  = x + T / B[I*l + i]  ;
            T  = x + T / B[L*i + l]  ;
            T  = max(T, minval);
            //PK = P[C*d + t] * K[I*d + i];
            PK = P[C*d + t] * K[D*i + d];
            f += PK / T ;
            g -= PK / (T * T) ;
    }}
    
    u  = - f * f / g ;
    v  = - f / g - x ;
    xp = u / c - v ;
    
    x = clamp(xp, minval, maxval) ;
}

b[L*d + l] = x;
}
