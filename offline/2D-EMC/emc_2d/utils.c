constant sampler_t trilinear = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_LINEAR ;
constant sampler_t nearest = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_NEAREST ;


// T[i, t, r]    = w[d] W[i, t, r] + b[d] B[i]
// logR[d, t, r] = beta sum_i K[i, d] log T[i, t, r] - T[i, t, r]

__kernel void calculate_LR (
    image2d_array_t I, 
    global float *LR,  
    global unsigned char *K, 
    global float *w,
    global float *b, 
    global float *B, 
    global float *C, 
    global float *R, 
    global float *rx, 
    global float *ry, 
    const float beta, 
    const float i0,
    const float dx,
    const int pixels, 
    const int d0)
{
    int Kindex   = get_global_id(2);
    int class    = get_global_id(1);
    int rotation = get_global_id(0);
    
    int frames    = get_global_size(2);
    int classes   = get_global_size(1);
    int rotations = get_global_size(0);
        
    int frame = d0 + Kindex;

    float R_l[4];
    float T;
    float logR = 0.;
    
    int i;
    
    for (i=0; i<4; i++) {
        R_l[i] = R[4*rotation + i];
    }

    float4 coord ;
    float4 W;

    coord.z = class ;

    for (i=0; i<pixels; i++) {
        coord.y = i0 + (R_l[0] * rx[i] + R_l[1] * ry[i]) / dx + 0.5;
        coord.x = i0 + (R_l[2] * rx[i] + R_l[3] * ry[i]) / dx + 0.5;
        
        W = read_imagef(I, trilinear, coord);
        
        T = w[frame] * C[i] * W.x + b[frame] * B[i];
        if (T > 0) 
            logR += K[pixels * Kindex + i] * log(T) - T ;
        //if ((Kindex == 0) && (class == 0) && (rotation == 0) && (i > 150) && (i < 200))
        //    printf("%d %e\n", i, T);
    }

    LR[rotations * classes * Kindex + rotations * class + rotation] = beta * logR;
}





// Wsums[t, r] = sum_i W[t, r, i] 
//             = sum_i C[i] I[i(t, r, i)]
__kernel void calculate_tomogram_sums (
    global float *Wsums,
    image2d_array_t I, 
    global float *C, 
    global float *R, 
    global float *rx, 
    global float *ry, 
    const float i0,
    const float dx,
    const int pixels)
{
    int class = get_global_id(0);
    int rotation = get_global_id(1);

    int classes = get_global_size(0);
    int rotations = get_global_size(1);

    float R_l[4];
    float T;
    
    int i;
    
    for (i=0; i<4; i++) {
        R_l[i] = R[4*rotation + i];
    }

    float4 coord ;
    float4 W;

    coord.z = class ;

    T = 0.;
    for (i=0; i<pixels; i++) {
        coord.y = i0 + (R_l[0] * rx[i] + R_l[1] * ry[i]) / dx + 0.5;
        coord.x = i0 + (R_l[2] * rx[i] + R_l[3] * ry[i]) / dx + 0.5;
        
        W = read_imagef(I, trilinear, coord);
        
        T += C[i] * W.x;
    }
    
    // may not be coallesed write (but not big so don't worry)
    Wsums[rotations * class + rotation] = T;
}


__kernel void test (
image2d_array_t I, 
global float *out,  
global float *LR,  
global unsigned char *K, 
global float *w,
global float *b, 
global float *B, 
global float *R, 
global float *rx, 
global float *ry, 
const float beta, 
const float i0,
const float dx,
const int pixels)
{
int frame = get_global_id(0);
int class = get_global_id(1);
int rotation = get_global_id(2);

int frames = get_global_size(0);
int classes = get_global_size(1);
int rotations = get_global_size(2);

float R_l[4];
float T;
float logR = 0.;

int i;

for (i=0; i<4; i++) {
    R_l[i] = R[4*rotation + i];
}

float4 coord ;
float4 W;

coord.z = class ;


for (i=0; i<pixels; i++) {
    coord.x = i0 + (R_l[0] * rx[i] + R_l[1] * ry[i]) / dx + 0.5;
    coord.y = i0 + (R_l[2] * rx[i] + R_l[3] * ry[i]) / dx + 0.5;
    
    W = read_imagef(I, trilinear, coord);
    
    out[pixels * rotations * class + pixels * rotation + i] = W.x ;
    
}

}


//loop over (t, r) with pixel-chunking 
//    - sparsity: find frames with low P value
//    load P[:, t, r] # small
//    c       = sum_d w[d] P[d]
//    xmax[i] = sum_d P[d] K[d, i] / c
//    
//    loop iters:
//        calculate W[i] <-- I
//        loop d :
//            T[i]    = W[i] + b[d] B[i] / w[d]
//            PK      = P[d] K[d, i] 
//            f[i]   += PK / T[i]
//            g[i]   += PK / T[i]^2
//        
//        step[i] = f[i] / g[i] * (1 - f[i] / c)
//            
//        W[i] += step[i]
            

__kernel void calculate_xmax_W (
global float *Wout,  
global unsigned char *K, 
global float *P, 
global int *frame_list,
const float c,
const int frames,
const int pixels)
{
int j, d;
float xmax = 0.;

int i = get_global_id(0);


for (j=0; j<frames; j++){
    d = frame_list[j];
    xmax += P[d] * K[d * pixels + i];
}
xmax /= c;

Wout[i] = xmax ;
}


__kernel void update_W (
global float *Wout,  
global float *B,  
global float *w,
global float *b, 
global unsigned char *K, 
global float *P, 
global int *frame_list,
const float c,
const int iters,
const int frames,
const int pixels)
{
int j, d, iter;
float xmax = 0.;
float T, f, g, step, PK;

int i = get_global_id(0);


for (j=0; j<frames; j++){
    d = frame_list[j];
    xmax += P[d] * K[d * pixels + i];
}
xmax /= c;

float W = xmax/2.;

for (iter=0; iter<iters; iter++){
    f = 0.;
    g = 0.;
    for (j=0; j<frames; j++){
        d = frame_list[j];
        T  = W + b[d] * B[i] / w[d];
        PK = P[d] * K[d * pixels + i];
        f += PK / T ;
        g -= PK / (T*T) ;
    }
    
    step = f / g * (1 - f / c);
    
    W += step;
    W = clamp(W, (float)1e-8, xmax) ;
}

Wout[i] = W;
}










//    each worker has a d-index
//    w = w[d]
//    b = b[d]
//    loop iters:
//        loop t,r,i :
//            W[t, r, i] <-- I, C 
//            T     = w + b B[i] / W[t, r, i]
//            PK    = P[t, r, d] K[i, d]
//            f[d] += PK / T
//            g[d] -= PK / T^2

kernel void update_w(
    image2d_array_t I, 
    global float* B, 
    global float* w, 
    global float* b, 
    global unsigned char* K, 
    global float* P, 
    global float* c, 
    global float* xmax, 
    const int   iters, 
    const int   frames, 
    const int   classes,
    const int   rotations,
    const int   pixels,
    const int   d,
    global float* C, 
    global float* R, 
    global float* rx, 
    global float* ry, 
    global int* class_rotation_list,
    const int class_rotation_list_len,
    const float i0,
    const float dx)
{
    //int d    = get_group_id(0);
    int wid  = get_local_id(0);
    int size = get_local_size(0);


    int i, t, r, j, iter;
    float c_l, b_l, xmax_l, T, step, PK;
    unsigned char K_l;

    local float f[256];
    local float g[256];
    
    local float x ;


    x      = w[d];
    b_l    = b[d];
    xmax_l = xmax[d];
    c_l    = c[d];

    float B_l;
    float fsum, gsum;

    float4 coord ;
    float4 W;
    float rx_l, ry_l, C_l;

    for (iter=0; iter<iters; iter++){
        f[wid] = 0.;
        g[wid] = 0.;
        for (i=wid; i<pixels; i+=size){
            B_l  = b_l * B[i] ;
            K_l  = K[i];
            rx_l = rx[i];
            ry_l = ry[i];
            C_l  = C[i];
        
        for (j = 0; j < class_rotation_list_len; j++){
            t = class_rotation_list[2*j + 0]; 
            r = class_rotation_list[2*j + 1]; 
              
            coord.y = i0 + (R[4*r + 0] * rx_l + R[4*r + 1] * ry_l) / dx + 0.5;
            coord.x = i0 + (R[4*r + 2] * rx_l + R[4*r + 3] * ry_l) / dx + 0.5;
            coord.z = t ;
            
            W = read_imagef(I, trilinear, coord);
            
            T = x + B_l / (C_l * W.x);
            
            PK      = P[t * rotations + r] * K_l;
            f[wid] += PK / T ;
            g[wid] -= PK / (T*T) ;
        }}
        
        // work group reduce
        barrier(CLK_LOCAL_MEM_FENCE);
        if (wid == 0) {
            fsum = 0 ;
            gsum = 0 ;
            for (i=0; i<size; i++) {
                fsum += f[i];
                gsum += g[i];
            }
        step = fsum / gsum * (1 - fsum / c_l);
        x   += step;
        x    = clamp(x, (float)1e-8, xmax_l) ;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        
    }

    if (wid == 0) w[d] = x;
}

__kernel void calculate_fg_w (
    image2d_array_t I, 
    global float *P,  
    global unsigned char *K, 
    global float *w,
    global float *b, 
    global float *B, 
    global float *C, 
    global float *R, 
    global float *rx, 
    global float *ry, 
    global float *f_g, 
    global float *g_g, 
    
    const float i0,
    const float dx,
    const int pixels, 
    const int d0)
{
    // frame index in sub-chunk
    int Kindex   = get_global_id(2);
    int class    = get_global_id(1);
    int rotation = get_global_id(0);
    
    int frames    = get_global_size(2);
    int classes   = get_global_size(1);
    int rotations = get_global_size(0);
        
    // frame index in MPI-chunk
    int frame = d0 + Kindex;
    
    float R_l[4];
    float T, PK;
    
    int i;
    
    for (i=0; i<4; i++) {
        R_l[i] = R[4*rotation + i];
    }
    
    float4 coord ;
    float4 W;
    
    coord.z = class ;
    
    float f = 0.;
    float g = 0.;
    for (i=0; i<pixels; i++) {
        coord.y = i0 + (R_l[0] * rx[i] + R_l[1] * ry[i]) / dx + 0.5;
        coord.x = i0 + (R_l[2] * rx[i] + R_l[3] * ry[i]) / dx + 0.5;
        
        W = read_imagef(I, trilinear, coord);
        
        // hopefully caching deals with repeated reads of w[frame], b[frame] and P over pixels
        T   = w[Kindex] + b[Kindex] * B[i]/ (C[i] * W.x) ;
        PK  = K[pixels * Kindex + i] * P[rotations * classes * frame + rotations * class + rotation] ;
        
        f += PK / T ;
        g -= PK / (T*T) ;
        
        
    }
    
    f_g[classes * rotations * Kindex + rotations * class + rotation ] = f;
    g_g[classes * rotations * Kindex + rotations * class + rotation ] = g;
}

__kernel void calculate_fg_b (
    image2d_array_t I, 
    global float *P,  
    global unsigned char *K, 
    global float *w,
    global float *b, 
    global float *B, 
    global float *C, 
    global float *R, 
    global float *rx, 
    global float *ry, 
    global float *f_g, 
    global float *g_g, 
    
    const float i0,
    const float dx,
    const int pixels, 
    const int d0)
{
    // frame index in sub-chunk
    int Kindex   = get_global_id(2);
    int class    = get_global_id(1);
    int rotation = get_global_id(0);
    
    int frames    = get_global_size(2);
    int classes   = get_global_size(1);
    int rotations = get_global_size(0);
        
    // frame index in MPI-chunk
    int frame = d0 + Kindex;
    
    float R_l[4];
    float T, PK;
    
    int i;
    
    for (i=0; i<4; i++) {
        R_l[i] = R[4*rotation + i];
    }
    
    float4 coord ;
    float4 W;
    
    coord.z = class ;
    
    float f = 0.;
    float g = 0.;
    for (i=0; i<pixels; i++) {
        coord.y = i0 + (R_l[0] * rx[i] + R_l[1] * ry[i]) / dx + 0.5;
        coord.x = i0 + (R_l[2] * rx[i] + R_l[3] * ry[i]) / dx + 0.5;
        
        W = read_imagef(I, trilinear, coord);
        
        // hopefully caching deals with repeated reads of w[frame], b[frame] and P over pixels
        T   = b[Kindex] + w[Kindex] * C[i] * W.x / B[i] ;
        PK  = K[pixels * Kindex + i] * P[rotations * classes * frame + rotations * class + rotation] ;
        
        f += PK / T ;
        g -= PK / (T*T) ;
    }
    
    f_g[classes * rotations * Kindex + rotations * class + rotation ] = f;
    g_g[classes * rotations * Kindex + rotations * class + rotation ] = g;
}




kernel void test_1d_image(
    image1d_t I, 
    global int* out)
{
    int i = get_global_id(0);
    
    int4 v ;

    v = read_imagei(I, i);
    out[i] = v.x;
}
