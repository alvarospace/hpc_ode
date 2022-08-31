/*
    -- MAGMA (version 2.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date April 2022

       @author Mark Gates
       @author Azzam Haidar
       
       @generated from src/zlarfb_gpu_gemm.cpp, normal z -> s, Wed Apr 20 17:36:17 2022
*/
#include "magma_internal.h"

/***************************************************************************//**
    Purpose
    -------
    SLARFB applies a real block reflector H or its transpose H^H to a
    REAL m by n matrix C, from the left.
    
    __Note that this function assumes__ that the upper part of dV is 0
    because it is referenced. Same for upper/lower part of dT.

    Arguments
    ---------
    @param[in]
    side    magma_side_t
      -     = MagmaLeft:      apply H or H^H from the Left
      -     = MagmaRight:     apply H or H^H from the Right

    @param[in]
    trans   magma_trans_t
      -     = MagmaNoTrans:    apply H   (No transpose)
      -     = MagmaTrans: apply H^H (Conjugate transpose)

    @param[in]
    direct  magma_direct_t
            Indicates how H is formed from a product of elementary
            reflectors
      -     = MagmaForward:  H = H(1) H(2) . . . H(k) (Forward)
      -     = MagmaBackward: H = H(k) . . . H(2) H(1) (Backward)

    @param[in]
    storev  magma_storev_t
            Indicates how the vectors which define the elementary
            reflectors are stored:
      -     = MagmaColumnwise: Columnwise
      -     = MagmaRowwise:    Rowwise

    @param[in]
    m       INTEGER
            The number of rows of the matrix C.

    @param[in]
    n       INTEGER
            The number of columns of the matrix C.

    @param[in]
    k       INTEGER
            The order of the matrix T (= the number of elementary
            reflectors whose product defines the block reflector).

    @param[in]
    dV      REAL array on the GPU, dimension
                (LDDV,K) if STOREV = MagmaColumnwise
                (LDDV,M) if STOREV = MagmaRowwise and SIDE = MagmaLeft
                (LDDV,N) if STOREV = MagmaRowwise and SIDE = MagmaRight
            The matrix V. See further details.

    @param[in]
    lddv    INTEGER
            The leading dimension of the array V.
            If STOREV = MagmaColumnwise and SIDE = MagmaLeft,  LDDV >= max(1,M);
            if STOREV = MagmaColumnwise and SIDE = MagmaRight, LDDV >= max(1,N);
            if STOREV = MagmaRowwise, LDDV >= K.

    @param[in]
    dT      REAL array on the GPU, dimension (LDDT,K)
            The triangular k by k matrix T in the representation of the
            block reflector.

    @param[in]
    lddt    INTEGER
            The leading dimension of the array T. LDDT >= K.

    @param[in,out]
    dC      REAL array on the GPU, dimension (LDDC,N)
            On entry, the m by n matrix C.
            On exit, C is overwritten by H*C, or H^H*C, or C*H, or C*H^H.

    @param[in]
    lddc    INTEGER
            The leading dimension of the array C. LDDC >= max(1,M).

    @param
    dwork   (workspace) REAL array, dimension (LDWORK,K)

    @param[in]
    ldwork  INTEGER
            The leading dimension of the array WORK.
            If SIDE = MagmaLeft,  LDWORK >= max(1,N);
            if SIDE = MagmaRight, LDWORK >= max(1,M);

    @param
    dworkvt (workspace) REAL array, dimension (LDWORKT,K)

    @param[in]
    ldworkvt INTEGER
            The leading dimension of the array WORKVT.
            LDWORKVT >= max(1,min(M,N));

    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    Further Details
    ---------------
    The shape of the matrix V and the storage of the vectors which define
    the H(i) is best illustrated by the following example with n = 5 and
    k = 3.
    All elements including 0's and 1's are stored, unlike LAPACK.

        DIRECT = MagmaForward and         DIRECT = MagmaForward and
        STOREV = MagmaColumnwise:         STOREV = MagmaRowwise:

                 V = (  1  0  0 )                 V = (  1 v1 v1 v1 v1 )
                     ( v1  1  0 )                     (  0  1 v2 v2 v2 )
                     ( v1 v2  1 )                     (  0  0  1 v3 v3 )
                     ( v1 v2 v3 )
                     ( v1 v2 v3 )

        DIRECT = MagmaBackward and        DIRECT = MagmaBackward and
        STOREV = MagmaColumnwise:         STOREV = MagmaRowwise:

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1  0  0 )
                     ( v1 v2 v3 )                     ( v2 v2 v2  1  0 )
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
                     (  0  1 v3 )
                     (  0  0  1 )

    @ingroup magma_larfb
*******************************************************************************/
extern "C" magma_int_t
magma_slarfb_gpu_gemm(
    magma_side_t side, magma_trans_t trans, magma_direct_t direct, magma_storev_t storev,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloat_const_ptr dV,    magma_int_t lddv,
    magmaFloat_const_ptr dT,    magma_int_t lddt,
    magmaFloat_ptr dC,          magma_int_t lddc,
    magmaFloat_ptr dwork,       magma_int_t ldwork,
    magmaFloat_ptr dworkvt,     magma_int_t ldworkvt,
    magma_queue_t queue )
{
    /* Constants */
    const float c_zero    = MAGMA_S_ZERO;
    const float c_one     = MAGMA_S_ONE;
    const float c_neg_one = MAGMA_S_NEG_ONE;

    magma_int_t info = 0;
    
    /* Function Body */
    if (m <= 0 || n <= 0) {
        return info;
    }
    
    /* Local variables */
    magma_int_t ldwvt = (m > n ?  k : m);
    magma_int_t ldw;
    if ( side == MagmaLeft ) {
        ldw = k;
    } else {
        ldw = m;
    }
    
    // whether V is stored transposed or not
    magma_trans_t notransV, transV;
    if (storev == MagmaColumnwise) {
        notransV = MagmaNoTrans;
        transV   = MagmaTrans;
    }
    else {
        notransV = MagmaTrans;
        transV   = MagmaNoTrans;
    }

    if ( side == MagmaLeft ) {
        // Form H C or H^H C
        // Comments assume H C.
        // When forming H^H C, T gets transposed via trans.
        
        // W = V^H C
        magma_sgemm( transV, MagmaNoTrans,
                     k, n, m,
                     c_one,  dV,    lddv,
                             dC,    lddc,
                     c_zero, dwork, ldw, queue );

        if (m <= n) {
            // W2 = V T
            magma_sgemm( notransV, trans,
                         m, k, k,
                         c_one,  dV, lddv,
                                 dT, lddt,
                         c_zero, dworkvt, ldwvt, queue );
            
            // C = C - W2 W = C - V T V^H C = (I - V T V^H) C = H C
            magma_sgemm( MagmaNoTrans, MagmaNoTrans,
                         m, n, k,
                         c_neg_one, dworkvt,  ldwvt,
                                    dwork,    ldw,
                         c_one,     dC,       lddc, queue );
        }
        else {
            // W2 = T W  = T V^H C
            magma_sgemm( trans, MagmaNoTrans,
                         k, n, k,
                         c_one,  dT,      lddt,
                                 dwork,   ldw,
                         c_zero, dworkvt, ldwvt, queue );
            
            // C = C - V W2 = C - V T V^H C = (I - V T V^H) C = H C
            magma_sgemm( notransV, MagmaNoTrans,
                         m, n, k,
                         c_neg_one, dV,      lddv,
                                    dworkvt, ldwvt,
                         c_one,     dC,      lddc, queue );
        }
    }
    else {
        // Form C H or C H^H
        // Comments assume C H.
        // When forming C H^H, T gets transposed via trans.
        
        // W = C V
        magma_sgemm( MagmaNoTrans, notransV,
                     m, k, n,
                     c_one,  dC,    lddc,
                             dV,    lddv,
                     c_zero, dwork, ldw, queue );
        
        if (m <= n) {
            // W2 = W T = C V T
            magma_sgemm( MagmaNoTrans, trans,
                         m, k, k,
                         c_one,  dwork,   ldw,
                                 dT,      lddt,
                         c_zero, dworkvt, ldwvt, queue );
            
            // C = C - W2 V^H = C - C V T V^H = C (I - V T V^H) = C H
            magma_sgemm( MagmaNoTrans, transV,
                         m, n, k,
                         c_neg_one, dworkvt, ldwvt,
                                    dV,      lddv,
                         c_one,     dC,      lddc, queue );
        }
        else {
            // W2 = T V^H
            magma_sgemm( trans, transV,
                         k, n, k,
                         c_one,  dT,      lddt,
                                 dV,      lddv,
                         c_zero, dworkvt, ldwvt, queue );
            
            // C = C - W W2 = C - C V T V^H = C (I - V T V^H) = C H
            magma_sgemm( MagmaNoTrans, MagmaNoTrans,
                         m, n, k,
                         c_neg_one, dwork,   ldw,
                                    dworkvt, ldwvt,
                         c_one,     dC,      lddc, queue );
        }
    }

    return info;
} /* magma_slarfb */
