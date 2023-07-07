//#include "stdio.h"

#define STENCIL(var) class Expr<double, void, opStencil>  var((double *)(_##var + i))

int sz2D, sz3D, szI, szJ, szK;

enum Op  {opValue, opStencil, opSum, opSub, opMul, opDiv, opRtoU, opRtoV,
          opDERtoU, opDNRtoV, opdivUVtoR};


extern "C" __global__ void initialize(unsigned int sizeK, unsigned int sizeJ, unsigned int sizeI)
{
    printf("%i %i %i \n", sizeK, sizeJ, sizeI);
    // Stores some global variables that are used a lot (instead of passing the as parameters on each function).
    szI  = sizeI;
    szJ  = sizeJ;
    szK  = sizeK;
    sz2D = szI*szJ;
    sz3D = szI*szJ*szK;
}


template<typename L, typename R, int op>
class Expr
{
};


template<typename L>
class Expr<L, void, opStencil>
{
    L * const p;


public:
    typedef L ResType;

    Expr(L * _p): p(_p)
    {
    }
    L &operator()(int const j, int const i) const
    {
        return *(p + j*szJ + i);
    }
    L operator=(L val) const
    {
//        printf("%p >-- %f\n",  p, expr.Eval(0,0));
        return (*p = val);
    }
    template<typename T1>
    void operator=(const T1 &expr) const
    {
////        *p = 1.0;
//        expr;
//        printf("%p -- %f\n",  p, 1.0); //*expr.Eval(0,0));
//expr.Eval(0,0);
////        printf("%p -sss- %f\n",  p, 1.0);

        *p = expr.Eval(0,0);
    }


//    void operator=(const double val) const
//    {
////        printf("%p ->- %f\n",  p, val); //expr.Eval(0,0));
//        *p = val;
//    }

    L Eval(int const j, int const i) const
    {
        return *(p + j*szJ + i);
    }




//    operator Expr<T, opStencil>() const { return Expr<T, opStencil>(*this); }
};




//template<typename T, typename L, typename R, int op>
//void operator=(Expr(T, void, opStencil) &l, Expr<L,R,op> const &expr) const
//{
//    printf("%p -- %f\n",  p, expr.Eval(0,0));
//    return (*l.p = expr.Eval(0,0));
//}




template<typename L>
class Expr<L, void, opRtoU>
{
    L &expr;

//    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;   // The product is a simple way to remove the &
public:
    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;

    Expr(L _expr): expr(_expr)
    {
    }

//    ~Expr()
//    {
//        printf("HERE!!!\n");
//    }

    ResType operator()(int const j, int const i) const
    {
        return this->Eval(i,j);
    }

    ResType Eval(int const j, int const i) const
    {
        return (expr.Eval(j, i) + expr.Eval(j-1, i)) * 0.5;
    }
};


template<typename L>
class Expr<L, void, opRtoV>
{
    L &expr;
//    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;   // The product is a simple way to remove the &
protected:
    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;
public:
    Expr(L _expr): expr(_expr)
    {
    }
    ResType operator()(int const j, int const i) const
    {
        return this->Eval(i,j);
    }

    ResType Eval(int const j, int const i) const
    {
        return (expr.Eval(j, i) + expr.Eval(j, i-1)) * 0.5;
    }

};


template<typename L>
class Expr<L, void, opDERtoU>
{
    L &expr;
//    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;   // The product is a simple way to remove the &

protected:
    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;
    typedef Expr<ResType, void, opStencil> OType;
    OType &on_u;

public:
    Expr(L _expr, OType _on_u): expr(_expr), on_u(_on_u)
    {
    }
    ResType operator()(int const j, int const i) const
    {
        return this->Eval(j,i);
    }

    ResType Eval(int const j, int const i) const
    {
        return on_u(j,i)*(expr.Eval(j,i) - expr.Eval(j,i-1));
    }

};


template<typename L>
class Expr<L, void, opDNRtoV>
{
    L &expr;
//    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;   // The product is a simple way to remove the &

protected:
    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;
    typedef Expr<ResType, void, opStencil> OType;
    OType &om_v;

public:
    Expr(L _expr, OType _om_v): expr(_expr), om_v(_om_v)
    {
    }
    ResType operator()(int const j, int const i) const
    {
        return this->Eval(j,i);
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
    typedef Expr<ResType, void, opStencil> OType;
    OType &on_u;
    OType &om_v;
    OType &pn;
    OType &pm;

public:
    Expr(L _U, R _V, OType _on_u, OType _om_v, OType _pn, OType _pm): U(_U), V(_V), on_u(_on_u), om_v(_om_v), pn(_pn), pm(_pm)
    {
    }
    ResType operator()(int const j, int const i) const
    {
        return this->Eval(j,i);
    }

    ResType Eval(int const j, int const i) const
    {
        return  ( (U.Eval(j, i+1)*on_u.Eval(j, i+1) - U.Eval(j, i)*on_u.Eval(j, i)) + (V.Eval(j+1, i)*om_v.Eval(j+1, i) - V.Eval(j, i)*om_v.Eval(j, i)) )*pm.Eval(j, i)*pn.Eval(j, i);
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
//        printf("-%f\n", val);
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
auto DERtoU(const T &R, const Expr<double, void, opStencil> &on_u)
{
    return Expr<T, void, opDERtoU>(R, on_u);
}

template<typename T>
auto DNRtoV(const T &R, const Expr<double, void, opStencil> &om_v)
{
    return Expr<T, void, opDNRtoV>(R, om_v);
}


template<typename T1, typename T2>
auto divUVtoR(const T1 &U, const T2 &V, const Expr<double, void, opStencil> &on_u, const Expr<double, void, opStencil> &om_v, const Expr<double, void, opStencil> &pn, const Expr<double, void, opStencil> &pm)
{
    return Expr<T1, T2, opdivUVtoR>(U, V, on_u, om_v, pn, pm);
}




extern "C"  __global__
void computeMomentumRHSPred(const double *_h, const double *_on_u, const double *_om_v,
                            const double *_rhs_ubar, const double *_rhs_vbar,
                            const double *_zeta_t1, const double *_zeta_t2, const double g)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    STENCIL(rhs_ubar);
    STENCIL(rhs_vbar);

//    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;
    if (((i % szI) == 0 || ((i % szI) == (szI - 1) || (i/szI) == 0) || (i/szI) == (szJ - 1)) || i >= sz2D)
    {
//        rhs_ubar = 0.0;
//        rhs_vbar = 0.0;
        return;
    }

    STENCIL(h);
    STENCIL(on_u);
    STENCIL(om_v);
    STENCIL(zeta_t1);
    STENCIL(zeta_t2);


    constexpr double weight = 2.0/5.0;
    auto gzeta  = (1 - weight)*zeta_t2 + weight*zeta_t1;

//    auto gzeta2 = gzeta*gzeta;
    rhs_ubar = 0.5*g*(RtoU(h)*DERtoU(gzeta,on_u) + DERtoU(gzeta,on_u));
    rhs_vbar = 0.5*g*(RtoV(h)*DNRtoV(gzeta,om_v) + DNRtoV(gzeta,om_v));
}


extern "C"  __global__
void computeMomentumRHSCorr(const double *_h, const double *_on_u, const double *_om_v,
                            const double *_rhs_ubar, const double *_rhs_vbar,
                            const double *_zeta_t0, const double *_zeta_t1, const double *_zeta_t2, const double g)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    STENCIL(rhs_ubar);
    STENCIL(rhs_vbar);

//    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;
    if (((i % szI) == 0 || ((i % szI) == (szI - 1) || (i/szI) == 0) || (i/szI) == (szJ - 1)) || i >= sz2D)
    {
//        rhs_ubar = 0.0;
//        rhs_vbar = 0.0;
        return;
    }

    STENCIL(h);
    STENCIL(on_u);
    STENCIL(om_v);
    STENCIL(zeta_t0);
    STENCIL(zeta_t1);
    STENCIL(zeta_t2);

    constexpr double weight = 4.0/25.0;
    auto gzeta = (1 - weight)*zeta_t1 + weight*0.5*(zeta_t2 + zeta_t0);

//    auto gzeta2 = gzeta*gzeta;
    rhs_ubar = 0.5*g*(RtoU(h)*DERtoU(gzeta,on_u) + DERtoU(gzeta,on_u));
    rhs_vbar = 0.5*g*(RtoV(h)*DNRtoV(gzeta,om_v) + DNRtoV(gzeta,om_v));
}


extern "C"  __global__
void computeZetaRHS(const double *_zeta, const double *_h, double *_ubar, const double *_vbar,
                    const double *_on_u, const double *_om_v, const double *_pn, const double *_pm, double *_res)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
//TODO: also (i % szJ) == N || (i/szJ) == M)
//    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;

    STENCIL(res);

    if (((i % szI) == 0 || ((i % szI) == (szI - 1) || (i/szI) == 0) || (i/szI) == (szJ - 1)) || i >= sz2D)
    {
//        res = 0.0;
        return;
    }

    STENCIL(zeta);
    STENCIL(h);
    STENCIL(ubar);
    STENCIL(vbar);
    STENCIL(on_u);
    STENCIL(om_v);
    STENCIL(pn);
    STENCIL(pm);


    // compute the water column depth
    auto D = zeta + h;

    auto DU = ubar*RtoU(D);
    auto DV = vbar*RtoV(D);

    res = divUVtoR(DU, DV, on_u, om_v, pn, pm);
}




extern "C"  __global__
void aaa(const double Dt, const double *_zeta_t0, const double *_zeta_t1, const double *_zeta_t2,
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
void bbb(const double *_zeta_t1, const double *_zeta_t2)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;

    STENCIL(zeta_t1);
    STENCIL(zeta_t2);


//    constexpr double weight = 2.0/5.0;
//    gzeta  = (1 - weight)*zeta_t2(0,0) + weight*zeta_t1(0,0);
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
void Pred2(const double Dt, const double *_u_t1, const double *_u_t2, const double *_v_t1, const double *_v_t2,
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