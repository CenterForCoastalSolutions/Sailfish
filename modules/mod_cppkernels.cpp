//#include "stdio.h"

#define STENCIL(var)   class Expr<double, void, opStencil2D>  var((double *)(_##var + i))
#define STENCIL3D(var) class Expr<double, void, opStencil3D>  var((double *)(_##var + i))


#define OPERATOR_ACCESS ResType operator()(int const k, int const j, int const i) const { return Eval(k,j,i);} \
                        ResType operator()(int const j, int const i) const { return Eval(j,i); }

int sz2D, sz3D, szI, szJ, szK;
double *_on_u, *_om_v, *_pm, *_pn;

enum Op  {opValue, opStencil2D, opStencil3D, opSum, opSub, opMul, opDiv, opRtoU, opRtoV,
          opDERtoU, opDXRtoV, opdivUVtoR,
          opDXXcentered, opUtoUW_4th};


extern "C" __global__ void initialize(unsigned int sizeK, unsigned int sizeJ, unsigned int sizeI,
                                      double *on_u, double *om_v, double *pn, double *pm)
{
    printf("Initialazing kernel information: %i %i %i \n", sizeK, sizeJ, sizeI);

    // Stores some global variables that are used a lot (instead of passing the as parameters on each function).
    szI  = sizeI;
    szJ  = sizeJ;
    szK  = sizeK;
    sz2D = szI*szJ;
    sz3D = szI*szJ*szK;

    _on_u = on_u;
    _om_v = om_v;
    _pn   = pn;
    _pm   = pm;

}


template<typename L, typename R, int op>
class Expr
{
};


template<typename L>
class Expr<L, void, opStencil2D>
{
    L * const p;


public:
    typedef L ResType;

    Expr(L * _p): p(_p)
    {
    }
    L &operator()(int const j, int const i) const
    {
        return *(p + j*szI + i);
    }

    // By design, we choose not to return a value from the assignment operator.
    void operator=(L val) const
    {
        *p = val;
    }

    // By design, we choose not to return a value from the assignment operator.
    template<typename T1>
    void operator=(const T1 &expr) const
    {
        *p = expr.Eval(0,0);
    }


    L Eval(int const j, int const i) const
    {
        return *(p + j*szI + i);
    }


};



template<typename L>
class Expr<L, void, opStencil3D>
{
    L * const p;


public:
    typedef L ResType;

    Expr(L * _p): p(_p)
    {
    }

    L &operator()(int const k, int const j, int const i) const
    {
        return Eval(k,j,i);
    }

    // By design, we choose not to return a value from the assignment operator.
    void operator=(L val) const
    {
        *p = val;
    }

    // By design, we choose not to return a value from the assignment operator.
    template<typename T1>
    void operator=(const T1 &expr) const
    {
        *p = expr.Eval(0,0,0);
    }

    L Eval(int const k, int const j, int const i) const
    {
        return *(p + k*sz2D + j*szI + i);
    }



};



template<typename L>
class Expr<L, void, opRtoU>
{
    L &expr;

public:
    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;

    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS

    ResType Eval(int const j, int const i) const
    {
        return (expr.Eval(j,i) + expr.Eval(j-1,i)) * 0.5;
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return (expr.Eval(k,j,i) + expr.Eval(k,j,i-1)) * 0.5;
    }
};


template<typename L>
class Expr<L, void, opRtoV>
{
    L &expr;
protected:
    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;
public:
    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS

    ResType Eval(int const j, int const i) const
    {
        return (expr.Eval(j, i) + expr.Eval(j-1, i)) * 0.5;
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return (expr.Eval(j, i) + expr.Eval(k,j-1, i)) * 0.5;
    }

};


template<typename L, typename R>
class Expr<L, R, opUtoUW_4th>
{
    L &Var;

protected:
    typedef decltype(Var.Eval(0,0)*Var.Eval(0,0)) ResType;

public:
    Expr(L _Var): Var(_Var)
    {
    }


    auto Eval(int const j, int const i) const { return (9.0/16.0)*(var(k,j,i) + var(k+1,j,i)) - (1.0/16.0)*(var(k-1,j,i) + var(k+2,j,i)); };
};



template<typename L>
class Expr<L, void, opDERtoU>
{
    L &expr;

protected:
    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;
    typedef Expr<ResType, void, opStencil2D> OType;
    OType &on_u;

public:
    Expr(L _expr, OType _on_u): expr(_expr), on_u(_on_u)
    {
    }

    ResType Eval(int const j, int const i) const
    {
        return on_u(j,i)*(expr.Eval(j,i) - expr.Eval(j,i-1));
    }

};


template<typename L>
class Expr<L, void, opDXRtoV>
{
    L &expr;

protected:
    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;
    typedef Expr<ResType, void, opStencil2D> OType;
    OType &om_v;

public:
    Expr(L _expr, OType _om_v): expr(_expr), om_v(_om_v)
    {
    }

    ResType Eval(int const j, int const i) const
    {
        return om_v(j,i)*(expr.Eval(j,i) - expr.Eval(j-1,i));
    }

};


template<typename L, typename R>
class Expr<L, R, opdivUVtoR>
{
    L &U;
    R &V;

protected:
    typedef decltype(U.Eval(0,0)*V.Eval(0,0)) ResType;
    typedef Expr<ResType, void, opStencil2D> OType;
    OType &on_u;
    OType &om_v;
    OType &pn;
    OType &pm;

public:
    Expr(L _U, R _V, OType _on_u, OType _om_v, OType _pn, OType _pm): U(_U), V(_V), on_u(_on_u), om_v(_om_v), pn(_pn), pm(_pm)
    {
    }

    ResType Eval(int const j, int const i) const
    {
        return  ( (U.Eval(j, i+1)*on_u.Eval(j, i+1) - U.Eval(j, i)*on_u.Eval(j, i)) + (V.Eval(j+1, i)*om_v.Eval(j+1, i) - V.Eval(j, i)*om_v.Eval(j, i)) )*pm.Eval(j, i)*pn.Eval(j, i);
    }

};


template<typename L>
class Expr<L, void, opDXXcentered>
{
    L &Var;

protected:
    typedef decltype(Var.Eval(0,0)*Var.Eval(0,0)) ResType;

public:
    Expr(L _Var): Var(_Var)
    {
    }

    ResType Eval(int const j, int const i) const
    {
        return Var.Eval(k,j,i-1) - 2.0*Var.Eval(k,j,i) + Var.Eval(k,j,i+1);
    }

};



template<typename L>
class Expr<L, void, opValue>
{
    const L  val;
public:
    Expr(const L _val): val(_val)  {}

    Expr(L &_val): val(_val)  {}

    const L Eval(int const j, int const i) const
    {
        return val;
    }


};




template<typename L, typename R>
class Expr<L, R, opSum>
{
    const L l;
    const R r;

    typedef decltype(l.Eval(0,0)*r.Eval(0,0)) ResType;

public:
    Expr(const L &_l, const R  &_r): l(_l), r(_r) {};

    ResType Eval(int const j, int const i) const { return l.Eval(j,i) + r.Eval(j,i); };
};


template<typename L, typename R>
class Expr<L, R, opSub>
{
    const L l;
    const R r;

    typedef decltype(l.Eval(0,0)*r.Eval(0,0)) ResType;

public:
    Expr(const L &_l, const R  &_r): l(_l), r(_r) {};

    ResType Eval(int const j, int const i) const { return l.Eval(j,i) - r.Eval(j,i); };
};


template<typename L, typename R>
class Expr<L, R, opMul>
{
    const L l;
    const R r;
    typedef decltype(l.Eval(0,0)*r.Eval(0,0)) ResType;

public:
    Expr(const L &_l, const R  &_r): l(_l), r(_r) {};

    auto Eval(int const j, int const i) const { return l.Eval(j,i) * r.Eval(j,i); };
};


template<typename L, typename R>
class Expr<L, R, opDiv>
{
    const L l;
    const R r;
    typedef decltype(l.Eval(0,0)*r.Eval(0,0)) ResType;

public:
    Expr(const L &_l, const R  &_r): l(_l), r(_r) {};

    auto Eval(int const j, int const i) const { return l.Eval(j,i) / r.Eval(j,i); };
};







template<typename L, typename R>
auto operator+(const L &l, const R &r) { return Expr<L, R, opSum>(l, r); }

template<typename L, typename R>
auto operator-(const L &l, const R &r) { return Expr<L, R, opSub>(l, r); }

template<typename L, typename R>
auto operator*(const L &l, const R &r) { return Expr<L, R, opMul>(l, r); }

template<typename L, typename R>
auto operator/(const L &l, const R &r) { return Expr<L, R, opDiv>(l, r); }




template<typename R>
auto operator+(const double l, const R &r) { auto ll = Expr<double, void, opValue>(l); return Expr<decltype(ll), R, opSum>(ll, r); }

template<typename R>
auto operator-(const double l, const R &r) { auto ll = Expr<double, void, opValue>(l); return Expr<decltype(ll), R, opSub>(ll, r); }

template<typename R>
auto operator*(const double l, const R &r) { auto ll = Expr<double, void, opValue>(l); return Expr<decltype(ll), R, opMul>(ll, r); }

template<typename R>
auto operator/(const double l, const R &r) { auto ll = Expr<double, void, opValue>(l); return Expr<decltype(ll), R, opDiv>(ll, r); }




template<typename L>
auto operator+(const L &l, const double r) { auto rr = Expr<double, void, opValue>(r); return Expr<L, decltype(rr), opSum>(l, rr); }

template<typename L>
auto operator-(const L &l, const double r) { auto rr = Expr<double, void, opValue>(r); return Expr<L, decltype(rr), opSub>(l, rr); }

template<typename L>
auto operator*(const L &l, const double r) { auto rr = Expr<double, void, opValue>(r); return Expr<L, decltype(rr), opMul>(l, rr); }

template<typename L>
auto operator/(const L &l, const double r) { auto rr = Expr<double, void, opValue>(r); return Expr<L, decltype(rr), opDiv>(l, rr); }


template<typename T>
auto RtoU(const T &R)
{
    return Expr<T, void, opRtoU>(R);
}


template<typename T>
auto RtoV(const T &R)
{
    return Expr<T, void, opRtoV>(R);
}


template<typename T>
auto UtoUW_4th(const T &V) {
    return Expr<T, void, opUtoUW_4th>(U);
}


template<typename T>
auto DERtoU(const T &R, const Expr<double, void, opStencil2D> &on_u)
{
    return Expr<T, void, opDERtoU>(R, on_u);
}

template<typename T>
auto DXRtoV(const T &R, const Expr<double, void, opStencil2D> &om_v)
{
    return Expr<T, void, opDXRtoV>(R, om_v);
}


template<typename T1, typename T2>
auto divUVtoR(const T1 &U, const T2 &V, const Expr<double, void, opStencil2D> &on_u, const Expr<double, void, opStencil2D> &om_v, const Expr<double, void, opStencil2D> &pn, const Expr<double, void, opStencil2D> &pm)
{
    return Expr<T1, T2, opdivUVtoR>(U, V, on_u, om_v, pn, pm);
}


template<typename T>
auto DXXUtoU(const T &U) { return Expr<T, void, opDXXcentered>(U); }

template<typename T>
auto DXXVtoV(const T &V) { return Expr<T, void, opDXXcentered>(U); }










extern "C"  __global__
void computeMomentumRHSPred(const double *_h,
                            const double *_rhs_ubar, const double *_rhs_vbar,
                            const double *_zeta_t1, const double *_zeta_t2, const double g, const double weight)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

//    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;
    if (((i % szI) == 0 || ((i % szI) == (szI - 1) || (i/szI) == 0) || (i/szI) == (szJ - 1)) || i >= sz2D)
    {
//        rhs_ubar = 0.0;
//        rhs_vbar = 0.0;
        return;
    }

    STENCIL(rhs_ubar);
    STENCIL(rhs_vbar);
    STENCIL(h);
    STENCIL(on_u);
    STENCIL(om_v);
    STENCIL(zeta_t1);
    STENCIL(zeta_t2);

    auto gzeta  = (1 - weight)*zeta_t2 + weight*zeta_t1;
    auto gzeta2 = gzeta*gzeta;   // TODO : sqr expression.

    rhs_ubar = 0.5*g*(RtoU(h)*DERtoU(gzeta,on_u) + DERtoU(gzeta2,on_u));
    rhs_vbar = 0.5*g*(RtoV(h)*DXRtoV(gzeta,om_v) + DXRtoV(gzeta2,om_v));

}


extern "C"  __global__
void computeMomentumRHSCorr(const double *_h,
                            const double *_rhs_ubar, const double *_rhs_vbar,
                            const double *_zeta_t0, const double *_zeta_t1, const double *_zeta_t2, const double g)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;



//    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;
    if (((i % szI) == 0 || ((i % szI) == (szI - 1) || (i/szI) == 0) || (i/szI) == (szJ - 1)) || i >= sz2D)
    {
//        rhs_ubar = 0.0;
//        rhs_vbar = 0.0;
        return;
    }

    STENCIL(rhs_ubar);
    STENCIL(rhs_vbar);
    STENCIL(h);
    STENCIL(on_u);
    STENCIL(om_v);
    STENCIL(zeta_t0);
    STENCIL(zeta_t1);
    STENCIL(zeta_t2);

    constexpr double weight = 4.0/25.0;
    auto gzeta = (1 - weight)*zeta_t1 + weight*0.5*(zeta_t2 + zeta_t0);

    auto gzeta2 = gzeta*gzeta;
    rhs_ubar = 0.5*g*(RtoU(h)*DERtoU(gzeta,on_u) + DERtoU(gzeta2,on_u));
    rhs_vbar = 0.5*g*(RtoV(h)*DXRtoV(gzeta,om_v) + DXRtoV(gzeta2,om_v));
}


extern "C"  __global__
void computeZetaRHS(const double *_zeta, const double *_h, double *_ubar, const double *_vbar, double *_res)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
//TODO: also (i % szJ) == N || (i/szJ) == M)
//    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;

    if (((i % szI) == 0 || ((i % szI) == (szI - 1) || (i/szI) == 0) || (i/szI) == (szJ - 1)) || i >= sz2D)
    {

        return;
    }

    STENCIL(res);
    STENCIL(zeta);
    STENCIL(h);
    STENCIL(ubar);
    STENCIL(vbar);
    STENCIL(on_u);
    STENCIL(om_v);
    STENCIL(pn);
    STENCIL(pm);


    // Water column depth
    auto D = zeta + h;

    // Fluxes.
    auto DU = ubar*RtoU(D);
    auto DV = vbar*RtoV(D);

    res = divUVtoR(DU, DV, on_u, om_v, pn, pm);

}




extern "C"  __global__
void computeZetaPred(const double Dt, const double *_zeta_t0, const double *_zeta_t1, const double *_zeta_t2,
                     const double *_rzeta_t1)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (((i % szI) == 0 || (i/szI) == 0) || i >= sz2D) return;

    STENCIL(zeta_t0);
    STENCIL(zeta_t1);
    STENCIL(zeta_t2);
    STENCIL(rzeta_t1);

    zeta_t2 = zeta_t0(0,0) + 2.0*Dt*rzeta_t1(0,0);


}



extern "C"  __global__
void Pred(const double Dt, const double *_v_t1, const double *_v_t2, const double *_rhs, const double *_D_t1, const double *_D_t2)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;

    STENCIL(v_t1);
    STENCIL(v_t2);
    STENCIL(D_t1);
    STENCIL(D_t2);
    STENCIL(rhs);

    v_t2 = (v_t1(0,0)*D_t1(0,0) + Dt*rhs(0,0))/D_t2(0,0);

}

extern "C"  __global__
void computeMomentumPred(const double Dt, const double *_u_t1, const double *_u_t2, const double *_v_t1, const double *_v_t2,
                         const double *_rhsu, const double *_rhsv, const double *_h,
                         const double *_zeta_t1, const double *_zeta_t2)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;

    STENCIL(u_t1);
    STENCIL(u_t2);
    STENCIL(v_t1);
    STENCIL(v_t2);
    STENCIL(zeta_t1);
    STENCIL(zeta_t2);
    STENCIL(h);
    STENCIL(rhsu);
    STENCIL(rhsv);

    auto D_t1 = zeta_t1 + h;
    auto D_t2 = zeta_t2 + h;
    auto D_t1U = RtoU(D_t1);
    auto D_t1V = RtoV(D_t1);
    auto D_t2U = RtoU(D_t2);
    auto D_t2V = RtoV(D_t2);

    u_t2 = (u_t1(0,0)*D_t1U(0,0) + Dt*rhsu(0,0))/D_t2U(0,0);
    v_t2 = (v_t1(0,0)*D_t1V(0,0) + Dt*rhsv(0,0))/D_t2V(0,0);

}



extern "C"  __global__
void AdamsMoultonCorr3rd(const double Dt, const double *_v_t2, const double *_rhs_t0, const double *_rhs_t1, const double *_rhs_t2)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    // Adams-Moulton 3rd order coefficients
    constexpr double AM3_2 =  5.0 / 12.0;
    constexpr double AM3_1 =  8.0 / 12.0;
    constexpr double AM3_0 = -1.0 / 12.0;

    if (i >= sz2D) return;

    STENCIL(rhs_t0);
    STENCIL(rhs_t1);
    STENCIL(rhs_t2);
    STENCIL(v_t2);

    v_t2 = v_t2(0,0) + Dt*(AM3_2*rhs_t2(0,0) + AM3_1*rhs_t1(0,0) + AM3_0*rhs_t0(0,0));
}



extern "C"  __global__
void AdamsMoultonCorr3rd2(const double Dt, const double *_u_t2, const double *_v_t2,
                          const double *_rhsu_t0, const double *_rhsu_t1, const double *_rhsu_t2,
                          const double *_rhsv_t0, const double *_rhsv_t1, const double *_rhsv_t2,
                          const double *_h, const double *_zeta_t1, const double *_zeta_t2)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    // Adams-Moulton 3rd order coefficients
    constexpr double AM3_2 =  5.0 / 12.0;
    constexpr double AM3_1 =  8.0 / 12.0;
    constexpr double AM3_0 = -1.0 / 12.0;

    if (i >= sz2D) return;

    STENCIL(rhsu_t0);
    STENCIL(rhsu_t1);
    STENCIL(rhsu_t2);
    STENCIL(rhsv_t0);
    STENCIL(rhsv_t1);
    STENCIL(rhsv_t2);
    STENCIL(u_t2);
    STENCIL(v_t2);
    STENCIL(zeta_t1);
    STENCIL(zeta_t2);
    STENCIL(h);

    auto D_t1 = zeta_t1 + h;
    auto D_t2 = zeta_t2 + h;
    auto D_t1U = RtoU(D_t1);
    auto D_t1V = RtoV(D_t1);
    auto D_t2U = RtoU(D_t2);
    auto D_t2V = RtoV(D_t2);

    u_t2 = (u_t2(0,0)*D_t1U + Dt*(AM3_2*rhsu_t2(0,0) + AM3_1*rhsu_t1(0,0) + AM3_0*rhsu_t0(0,0)))/D_t2U;
    v_t2 = (v_t2(0,0)*D_t1V + Dt*(AM3_2*rhsv_t2(0,0) + AM3_1*rhsv_t1(0,0) + AM3_0*rhsv_t0(0,0)))/D_t2V;
}