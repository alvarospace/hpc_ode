/*
    -- MAGMA (version 2.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date April 2022

       @generated from testing/testing_ztrsm_batched.cpp, normal z -> c, Wed Apr 20 17:37:23 2022
       @author Chongxiao Cao
       @author Tingxing Dong
       @author Azzam Haidar
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "testings.h"

#include "../control/batched_kernel_param.h"  // internal header; for TRI_NB

#if defined(_OPENMP)
#include <omp.h>
#include "../control/magma_threadsetting.h"  // internal header
#endif

#define h_A(i,j,s) (h_A + (i) + (j)*lda + (s)*lda*Ak)


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing ctrsm_batched
*/
int main( int argc, char** argv)
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();

    real_Double_t   gflops, magma_perf, magma_time=0, cublas_perf=0, cublas_time=0, cpu_perf=0, cpu_time=0;
    float          error, magma_error, cublas_error, lapack_error, work[1];
    magma_int_t M, N, info;
    magma_int_t Ak;
    magma_int_t sizeA, sizeB;
    magma_int_t lda, ldb, ldda, lddb;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t *ipiv;

    magmaFloatComplex c_zero = MAGMA_C_ZERO;

    magmaFloatComplex *h_A, *h_B, *h_Bcublas, *h_Bmagma, *h_Blapack, *h_X;
    magmaFloatComplex *d_A, *d_B;
    magmaFloatComplex **h_A_array = NULL;
    magmaFloatComplex **h_B_array = NULL;
    magmaFloatComplex **d_A_array = NULL;
    magmaFloatComplex **d_B_array = NULL;

    magmaFloatComplex **dW1_displ  = NULL;
    magmaFloatComplex **dW2_displ  = NULL;
    magmaFloatComplex **dW3_displ  = NULL;
    magmaFloatComplex **dW4_displ  = NULL;
    magmaFloatComplex **dinvA_array = NULL;
    magmaFloatComplex **dwork_array = NULL;

    magmaFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    magmaFloatComplex c_one = MAGMA_C_ONE;
    magmaFloatComplex alpha = MAGMA_C_MAKE(  0.29, -0.86 );
    int status = 0;
    magma_int_t batchCount;

    magma_opts opts( MagmaOptsBatched );
    opts.parse_opts( argc, argv );
    batchCount = opts.batchcount;

    float tol = opts.tolerance * lapackf77_slamch("E");
    printf("%% side = %s, uplo = %s, transA = %s, diag = %s \n",
           lapack_side_const(opts.side), lapack_uplo_const(opts.uplo),
           lapack_trans_const(opts.transA), lapack_diag_const(opts.diag) );
    printf("%% BatchCount   M     N   MAGMA Gflop/s (ms)   %s Gflop/s (ms)    CPU Gflop/s (ms)      MAGMA     %s   LAPACK error\n", g_platform_str, g_platform_str);
    printf("%%=========================================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            gflops = FLOPS_CTRSM(opts.side, M, N) / 1e9 * batchCount;

            if ( opts.side == MagmaLeft ) {
                lda = M;
                Ak  = M;
            } else {
                lda = N;
                Ak  = N;
            }

            ldb = M;

            ldda = magma_roundup( lda, opts.align );  // multiple of 32 by default
            lddb = magma_roundup( ldb, opts.align );  // multiple of 32 by default

            sizeA = lda*Ak*batchCount;
            sizeB = ldb*N*batchCount;
            magma_int_t NN = ldb*N;

            TESTING_CHECK( magma_cmalloc_cpu( &h_A,       sizeA  ));
            TESTING_CHECK( magma_cmalloc_cpu( &h_B,       sizeB   ));
            TESTING_CHECK( magma_cmalloc_cpu( &h_X,       sizeB   ));
            TESTING_CHECK( magma_cmalloc_cpu( &h_Blapack, sizeB   ));
            TESTING_CHECK( magma_cmalloc_cpu( &h_Bcublas, sizeB   ));
            TESTING_CHECK( magma_cmalloc_cpu( &h_Bmagma,  sizeB   ));
            TESTING_CHECK( magma_imalloc_cpu( &ipiv,      Ak      ));

            TESTING_CHECK( magma_cmalloc( &d_A,       ldda*Ak*batchCount ));
            TESTING_CHECK( magma_cmalloc( &d_B,       lddb*N*batchCount  ));

            TESTING_CHECK( magma_malloc_cpu( (void**) &h_A_array,   batchCount * sizeof(magmaFloatComplex*) ));
            TESTING_CHECK( magma_malloc_cpu( (void**) &h_B_array,   batchCount * sizeof(magmaFloatComplex*) ));

            TESTING_CHECK( magma_malloc( (void**) &d_A_array,   batchCount * sizeof(magmaFloatComplex*) ));
            TESTING_CHECK( magma_malloc( (void**) &d_B_array,   batchCount * sizeof(magmaFloatComplex*) ));

            magmaFloatComplex* dinvA=NULL;
            magmaFloatComplex* dwork=NULL; // invA and work are workspace in ctrsm
            magma_int_t dinvA_batchSize = magma_roundup( Ak, TRI_NB )*TRI_NB;
            magma_int_t dwork_batchSize = lddb*N;
            if(opts.version == 1) {
                // these allocations are needed only for the inversion-based trsm
                TESTING_CHECK( magma_malloc( (void**) &dW1_displ,   batchCount * sizeof(magmaFloatComplex*) ));
                TESTING_CHECK( magma_malloc( (void**) &dW2_displ,   batchCount * sizeof(magmaFloatComplex*) ));
                TESTING_CHECK( magma_malloc( (void**) &dW3_displ,   batchCount * sizeof(magmaFloatComplex*) ));
                TESTING_CHECK( magma_malloc( (void**) &dW4_displ,   batchCount * sizeof(magmaFloatComplex*) ));

                TESTING_CHECK( magma_malloc( (void**) &dinvA_array, batchCount * sizeof(magmaFloatComplex*) ));
                TESTING_CHECK( magma_malloc( (void**) &dwork_array, batchCount * sizeof(magmaFloatComplex*) ));
                TESTING_CHECK( magma_cmalloc( &dinvA, dinvA_batchSize * batchCount ));
                TESTING_CHECK( magma_cmalloc( &dwork, dwork_batchSize * batchCount ));

                magma_cset_pointer( dwork_array, dwork, lddb, 0, 0, dwork_batchSize, batchCount, opts.queue );
                magma_cset_pointer( dinvA_array, dinvA, magma_roundup( Ak, TRI_NB ), 0, 0, dinvA_batchSize, batchCount, opts.queue );
            }

            memset( h_Bmagma, 0, batchCount*ldb*N*sizeof(magmaFloatComplex) );
            magmablas_claset( MagmaFull, lddb, N*batchCount, c_zero, c_zero, dwork, lddb, opts.queue );

            /* Initialize the matrices */
            /* Factor A into LU to get well-conditioned triangular matrix.
             * Copy L to U, since L seems okay when used with non-unit diagonal
             * (i.e., from U), while U fails when used with unit diagonal. */
            lapackf77_clarnv( &ione, ISEED, &sizeA, h_A );

            for (int s=0; s < batchCount; s++) {
                lapackf77_cgetrf( &Ak, &Ak, h_A + s*lda*Ak, &lda, ipiv, &info );
                for( int j = 0; j < Ak; ++j ) {
                    for( int i = 0; i < j; ++i ) {
                        *h_A(i,j,s) = *h_A(j,i,s);
                    }
                }
            }

            lapackf77_clarnv( &ione, ISEED, &sizeB, h_B );
            memcpy( h_Blapack, h_B, sizeB*sizeof(magmaFloatComplex) );

            /* =====================================================================
               Performs operation using MAGMABLAS
               =================================================================== */
            magma_csetmatrix( Ak, Ak*batchCount, h_A, lda, d_A, ldda, opts.queue );
            magma_csetmatrix( M,  N*batchCount,  h_B, ldb, d_B, lddb, opts.queue );

            magma_cset_pointer( d_A_array, d_A, ldda, 0, 0, ldda*Ak, batchCount, opts.queue );
            magma_cset_pointer( d_B_array, d_B, lddb, 0, 0, lddb*N, batchCount, opts.queue );

            magma_time = magma_sync_wtime( opts.queue );
            if (opts.version == 1) {
                magmablas_ctrsm_inv_outofplace_batched(
                    opts.side, opts.uplo, opts.transA, opts.diag, 1,
                    M, N, alpha,
                    d_A_array,    ldda, // dA
                    d_B_array,    lddb, // dB
                    dwork_array,  lddb, // dX output
                    dinvA_array,  dinvA_batchSize,
                    dW1_displ,   dW2_displ,
                    dW3_displ,   dW4_displ,
                    1, batchCount, opts.queue);
                magma_time = magma_sync_wtime( opts.queue ) - magma_time;
                magma_perf = gflops / magma_time;
                magma_cgetmatrix( M, N*batchCount, dwork, lddb, h_Bmagma, ldb, opts.queue );
            }
            else {
                magmablas_ctrsm_batched(
                    opts.side, opts.uplo, opts.transA, opts.diag,
                    M, N, alpha,
                    d_A_array, ldda,
                    d_B_array, lddb,
                    batchCount, opts.queue );
                magma_time = magma_sync_wtime( opts.queue ) - magma_time;
                magma_perf = gflops / magma_time;
                magma_cgetmatrix( M, N*batchCount, d_B, lddb, h_Bmagma, ldb, opts.queue );
            }

            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_csetmatrix( M, N*batchCount, h_B, ldb, d_B, lddb, opts.queue );
            magma_cset_pointer( d_B_array, d_B, lddb, 0, 0, lddb*N, batchCount, opts.queue );

            cublas_time = magma_sync_wtime( opts.queue );
            #ifdef MAGMA_HAVE_CUDA
            cublasCtrsmBatched(
                opts.handle, cublas_side_const(opts.side), cublas_uplo_const(opts.uplo),
                cublas_trans_const(opts.transA), cublas_diag_const(opts.diag),
                int(M), int(N), (const cuFloatComplex*)&alpha,
                (const cuFloatComplex**) d_A_array, int(ldda),
                (      cuFloatComplex**) d_B_array, int(lddb), int(batchCount) );
            #else
            hipblasCtrsmBatched(
                opts.handle, cublas_side_const(opts.side), cublas_uplo_const(opts.uplo),
                cublas_trans_const(opts.transA), cublas_diag_const(opts.diag),
                int(M), int(N), (const hipblasComplex*)&alpha,
                (hipblasComplex* const*) d_A_array, int(ldda),
                (      hipblasComplex**) d_B_array, int(lddb), int(batchCount) );
            #endif

            cublas_time = magma_sync_wtime( opts.queue ) - cublas_time;
            cublas_perf = gflops / cublas_time;

            magma_cgetmatrix( M, N*batchCount, d_B, lddb, h_Bcublas, ldb, opts.queue );

            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                // populate pointer arrays on the host
                for(int s = 0; s < batchCount; s++){
                    h_A_array[s] = h_A + s * lda * Ak;
                    h_B_array[s] = h_Blapack + s * ldb * N;
                }
                cpu_time = magma_wtime();
                blas_ctrsm_batched(
                    opts.side, opts.uplo, opts.transA, opts.diag,
                    M, N,
                    alpha, h_A_array, lda,
                           h_B_array, ldb, batchCount );
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }

            /* =====================================================================
               Check the result
               =================================================================== */
            // ||b - 1/alpha*A*x|| / (||A||*||x||)
            magmaFloatComplex inv_alpha = MAGMA_C_DIV( c_one, alpha );
            float normR, normX, normA;
            magma_error  = 0;
            cublas_error = 0;

            memcpy( h_X, h_Bmagma, sizeB*sizeof(magmaFloatComplex) );

            // check magma
            for (int s=0; s < batchCount; s++) {
                normA = lapackf77_clantr( "M",
                                          lapack_uplo_const(opts.uplo),
                                          lapack_diag_const(opts.diag),
                                          &Ak, &Ak, h_A + s*lda*Ak, &lda, work );
                blasf77_ctrmm(
                    lapack_side_const(opts.side), lapack_uplo_const(opts.uplo),
                    lapack_trans_const(opts.transA), lapack_diag_const(opts.diag),
                    &M, &N, &inv_alpha,
                    h_A + s*lda*Ak, &lda,
                    h_X + s*ldb*N,  &ldb );

                blasf77_caxpy( &NN, &c_neg_one, h_B + s*ldb*N, &ione, h_X + s*ldb*N, &ione );

                normR = lapackf77_clange( "M", &M, &N, h_X      + s*ldb*N, &ldb, work );
                normX = lapackf77_clange( "M", &M, &N, h_Bmagma + s*ldb*N, &ldb, work );
                error = normR/(normX*normA);
                magma_error = magma_max_nan( error, magma_error );
            }

            memcpy( h_X, h_Bcublas, sizeB*sizeof(magmaFloatComplex) );
            // check cublas
            for (int s=0; s < batchCount; s++) {
                normA = lapackf77_clantr( "M",
                                          lapack_uplo_const(opts.uplo),
                                          lapack_diag_const(opts.diag),
                                          &Ak, &Ak, h_A + s*lda*Ak, &lda, work );
                blasf77_ctrmm(
                    lapack_side_const(opts.side), lapack_uplo_const(opts.uplo),
                    lapack_trans_const(opts.transA), lapack_diag_const(opts.diag),
                    &M, &N, &inv_alpha,
                    h_A + s*lda*Ak, &lda,
                    h_X + s*ldb*N, &ldb );

                blasf77_caxpy( &NN, &c_neg_one, h_B + s*ldb*N, &ione, h_X  + s*ldb*N, &ione );
                normR = lapackf77_clange( "M", &M, &N, h_X  + s*ldb*N,       &ldb, work );
                normX = lapackf77_clange( "M", &M, &N, h_Bcublas  + s*ldb*N, &ldb, work );
                error = normR/(normX*normA);
                cublas_error = magma_max_nan( error, cublas_error );
            }
            bool okay = (magma_error < tol && cublas_error < tol);
            status += ! okay;

            if ( opts.lapack ) {
                // check lapack
                // this verifies that the matrix wasn't so bad that it couldn't be solved accurately.
                lapack_error = 0.0;
                memcpy( h_X, h_Blapack, sizeB*sizeof(magmaFloatComplex) );
                for (int s=0; s < batchCount; s++) {
                    normA = lapackf77_clantr( "M",
                                              lapack_uplo_const(opts.uplo),
                                              lapack_diag_const(opts.diag),
                                              &Ak, &Ak, h_A + s*lda*Ak, &lda, work );
                    blasf77_ctrmm(
                        lapack_side_const(opts.side), lapack_uplo_const(opts.uplo),
                        lapack_trans_const(opts.transA), lapack_diag_const(opts.diag),
                        &M, &N, &inv_alpha,
                        h_A + s*lda*Ak, &lda,
                        h_X + s*ldb*N,  &ldb );

                    blasf77_caxpy( &NN, &c_neg_one, h_B + s*ldb*N, &ione, h_X + s*ldb*N, &ione );
                    normR = lapackf77_clange( "M", &M, &N, h_X + s*ldb*N,       &ldb, work );
                    normX = lapackf77_clange( "M", &M, &N, h_Blapack + s*ldb*N, &ldb, work );
                    error = normR/(normX*normA);
                    lapack_error = magma_max_nan( error, lapack_error );
                }

                printf("%10lld %5lld %5lld    %7.2f (%7.2f)     %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %8.2e   %8.2e   %s\n",
                        (long long) batchCount, (long long) M, (long long) N,
                        magma_perf,  1000.*magma_time,
                        cublas_perf, 1000.*cublas_time,
                        cpu_perf,    1000.*cpu_time,
                        magma_error, cublas_error, lapack_error,
                        (okay ? "ok" : "failed"));
            }
            else {
                printf("%10lld %5lld %5lld    %7.2f (%7.2f)     %7.2f (%7.2f)     ---   (  ---  )   %8.2e   %8.2e     ---      %s\n",
                        (long long) batchCount, (long long) M, (long long) N,
                        magma_perf,  1000.*magma_time,
                        cublas_perf, 1000.*cublas_time,
                        magma_error, cublas_error,
                        (okay ? "ok" : "failed"));
            }

            magma_free_cpu( h_A );
            magma_free_cpu( h_B );
            magma_free_cpu( h_X );
            magma_free_cpu( h_Blapack );
            magma_free_cpu( h_Bcublas );
            magma_free_cpu( h_Bmagma  );
            magma_free_cpu( ipiv );
            magma_free_cpu( h_A_array );
            magma_free_cpu( h_B_array );

            magma_free( d_A );
            magma_free( d_B );
            magma_free( d_A_array );
            magma_free( d_B_array );

            if(opts.version == 1) {
                magma_free( dW1_displ );
                magma_free( dW2_displ );
                magma_free( dW3_displ );
                magma_free( dW4_displ );

                magma_free( dinvA );
                magma_free( dwork );
                magma_free( dwork_array );
                magma_free( dinvA_array );
            }

            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    opts.cleanup();
    TESTING_CHECK( magma_finalize() );
    return status;
}
