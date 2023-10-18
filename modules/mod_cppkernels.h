

#define STENCIL(var)   class Expr<double, void, opStencil2D>  var((double *)(_##var + i))
#define STENCIL3D(var) class Expr<double, void, opStencil3D>  var((double *)(_##var + i))


#define OPERATOR_ACCESS ResType operator()(int const k, int const j, int const i) const { return Eval(k,j,i);} \
                        ResType operator()(int const j, int const i) const { return Eval(j,i); }

int sz2D, sz3D, szI, szJ, szK, K;
double *_on_u, *_om_v, *_pm, *_pn;

enum Op  {opValue, opStencil2D, opStencil3D, opSum, opSub, opMul, opDiv, opRtoU, opRtoV,
          opUtoR, opVtoR, opUtoP, opVtoP, opUpwindUtoR, opUpwindVtoR, opUpwindUtoP, opUpwindVtoP, opUtoUW_4th,
          opDxRtoU, opDyRtoV, opdivUVtoR,
          opDXXcentered, opDEEcentered, opDX, opDE};

enum nodeType {ntR, ntU, ntV, ntP};



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

    K = 0;

}

extern "C" __global__ void setK(unsigned int k)
{
    K = k;
}


template<class T>
bool zeroOf(T const &t) { return t.zero(); };


template<>
bool zeroOf(double const &t) { return 0.0; };




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

    template<typename T1>
    void operator-=(const T1 &expr) const
    {
        *p = *p - expr.Eval(0,0,0);
    }

    template<typename T1>
    void operator+=(const T1 &expr) const
    {
        *p = *p + expr.Eval(0,0,0);
    }


    L Eval(int const j, int const i) const
    {
        return *(p + j*szI + i);
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;

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
        *(p + K*sz2D) = val;
    }


    // By design, we choose not to return a value from the assignment operator.
    template<typename T1>
    void operator=(const T1 &expr) const
    {
        *(p + K*sz2D) = expr.Eval(0,0,0);
    }

    template<typename T1>
    void operator-=(const T1 &expr) const
    {
        *(p + K*sz2D) = *(p + K*sz2D) - expr.Eval(0,0,0);
    }

    template<typename T1>
    void operator+=(const T1 &expr) const
    {
        *(p + K*sz2D) = *(p + K*sz2D) + expr.Eval(0,0,0);
    }

    L Eval(int const k, int const j, int const i) const
    {
        return *(p + (k + K)*sz2D + j*szI + i);
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;

};


template<typename L>
class Expr<L, void, opUtoR>
{
    L &expr;

public:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;


    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS

    ResType Eval(int const j, int const i) const
    {
        return (expr.Eval(j,i) + expr.Eval(j,i+1))*0.5;
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return (expr.Eval(k,j,i) + expr.Eval(k,j,i+1))*0.5;
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L>
class Expr<L, void, opVtoR>
{
    L &expr;

public:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;


    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS

    ResType Eval(int const j, int const i) const
    {
        return (expr.Eval(j,i) + expr.Eval(j+1,i))*0.5;
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return (expr.Eval(k,j,i) + expr.Eval(k,j+1,i))*0.5;
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L>
class Expr<L, void, opUtoP>
{
    L &expr;

public:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;


    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS

    ResType Eval(int const j, int const i) const
    {
        return (expr.Eval(j,i) + expr.Eval(j-1,i))*0.5;
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return (expr.Eval(k,j,i) + expr.Eval(k,j-1,i))*0.5;
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L>
class Expr<L, void, opVtoP>
{
    L &expr;

public:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;


    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS

    ResType Eval(int const j, int const i) const
    {
        return (expr.Eval(j,i) + expr.Eval(j-1,i))*0.5;
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return (expr.Eval(k,j,i) + expr.Eval(k,j-1,i))*0.5;
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};



template<typename L>
class Expr<L, void, opRtoU>
{
    L &expr;

public:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS

    ResType Eval(int const j, int const i) const
    {
        return (expr.Eval(j,i) + expr.Eval(j,i-1))*0.5;
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return (expr.Eval(k,j,i) + expr.Eval(k,j,i-1))*0.5;
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L>
class Expr<L, void, opRtoV>
{
    L &expr;
protected:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

public:
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
        return (expr.Eval(k,j,i) + expr.Eval(k,j-1,i)) * 0.5;
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;

};


template<typename L, typename R>
class Expr<L, R, opUpwindUtoR>
{
    L &expr;
    R &vel;

protected:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

public:
    Expr(L _expr, R _vel): expr(_expr), vel(_vel)
    {
    }

    ResType Eval(int const k, int const j, int const i) const { (vel.Eval(k,j,i) + vel.Eval(k,j,i+1) >= 0)?expr.Eval(k,j,i):expr.Eval(k,j,i+1); };
    ResType Eval(int const j, int const i) const              { (vel.Eval(j,i)   + vel.Eval(j,i+1)   >= 0)?expr.Eval(j,i)  :expr.Eval(j,i+1); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};

template<typename L, typename R>
class Expr<L, R, opUpwindVtoR>
{
    L &expr;
    R &vel;

protected:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

public:
    Expr(L _expr, R _vel): expr(_expr), vel(_vel)
    {
    }

    ResType Eval(int const k, int const j, int const i) const { (vel.Eval(k,j,i) + vel.Eval(k,j+1,i) >= 0)?expr.Eval(k,j,i):expr.Eval(k,j+1,i); };
    ResType Eval(int const j, int const i) const              { (vel.Eval(j,i)   + vel.Eval(j+1,i)   >= 0)?expr.Eval(j,i)  :expr.Eval(j+1,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, typename R>
class Expr<L, R, opUpwindUtoP>
{
    L &expr;
    R &vel;

protected:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

public:
    Expr(L _expr, R _vel): expr(_expr), vel(_vel)
    {
    }

    ResType Eval(int const k, int const j, int const i) const { (vel.Eval(k,j,i) + vel.Eval(k,j-1,i) >= 0)?expr.Eval(k,j,i):expr.Eval(k,j-1,i); };
    ResType Eval(int const j, int const i) const              { (vel.Eval(j,i)   + vel.Eval(j-1,i)   >= 0)?expr.Eval(j,i)  :expr.Eval(j-1,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};

template<typename L, typename R>
class Expr<L, R, opUpwindVtoP>
{
    L &expr;
    R &vel;

protected:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

public:
    Expr(L _expr, R _vel): expr(_expr), vel(_vel)
    {
    }

    ResType Eval(int const k, int const j, int const i) const { (vel.Eval(k,j,i) + vel.Eval(k,j,i-1) >= 0)?expr.Eval(k,j,i):expr.Eval(k,j,i-1); };
    ResType Eval(int const j, int const i) const              { (vel.Eval(j,i)   + vel.Eval(j,i-1)   >= 0)?expr.Eval(j,i)  :expr.Eval(j,i-1); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};




template<typename L, typename R>
class Expr<L, R, opUtoUW_4th>
{
    L &Var;

protected:
    typedef decltype(zeroOf(Var)*zeroOf(Var)) ResType;

public:
    Expr(L _Var): Var(_Var)
    {
    }


    ResType Eval(int const k, int const j, int const i) const { return (9.0/16.0)*(var(k,j,i) + var(k+1,j,i)) - (1.0/16.0)*(var(k-1,j,i) + var(k+2,j,i)); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};






template<typename L>
class Expr<L, void, opDxRtoU>
{
    L &expr;

protected:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;
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

    ResType Eval(int const k, int const j, int const i) const
    {
        return on_u(j,i)*(expr.Eval(k,j,i) - expr.Eval(k,j,i-1));
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L>
class Expr<L, void, opDyRtoV>
{
    L &expr;

protected:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;
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

    ResType Eval(int const k, int const j, int const i) const
    {
        return om_v(j,i)*(expr.Eval(k,j,i) - expr.Eval(k,j-1,i));
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};



template<typename L>
class Expr<L, void, opDX>
{
    L &expr;

protected:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;
    typedef Expr<ResType, void, opStencil2D> OType;

public:
    Expr(L _expr): expr(_expr)
    {
    }

    ResType Eval(int const j, int const i) const
    {
        return expr.Eval(j,i) - expr.Eval(j,i-1);
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return expr.Eval(k,j,i) - expr.Eval(k,j,i-1);
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L>
class Expr<L, void, opDE>
{
    L &expr;

protected:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;
    typedef Expr<ResType, void, opStencil2D> OType;

public:
    Expr(L _expr): expr(_expr)
    {
    }

    ResType Eval(int const j, int const i) const
    {
        return expr.Eval(j,i) - expr.Eval(j-1,i);
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return expr.Eval(k,j,i) - expr.Eval(k,j-1,i);
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};



template<typename L, typename R>
class Expr<L, R, opdivUVtoR>
{
    L &U;
    R &V;

protected:
    typedef decltype(zeroOf(U)*zeroOf(V)) ResType;
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

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L>
class Expr<L, void, opDXXcentered>
{
    L &Var;

protected:
    typedef decltype(zeroOf(Var)*zeroOf(Var)) ResType;

public:
    Expr(L _Var): Var(_Var)
    {
    }

    ResType Eval(int const j, int const i) const
    {
        return Var.Eval(j,i-1) - 2.0*Var.Eval(j,i) + Var.Eval(j,i+1);
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return Var.Eval(k,j,i-1) - 2.0*Var.Eval(k,j,i) + Var.Eval(k,j,i+1);
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L>
class Expr<L, void, opDEEcentered>
{
    L &Var;

protected:
    typedef decltype(zeroOf(Var)*zeroOf(Var)) ResType;

public:
    Expr(L _Var): Var(_Var)
    {
    }

    ResType Eval(int const j, int const i) const
    {
        return Var.Eval(j-1,i) - 2.0*Var.Eval(j,i) + Var.Eval(j+1,i);
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return Var.Eval(k,j-1,i) - 2.0*Var.Eval(k,j,i) + Var.Eval(k,j+1,i);
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
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

    const L Eval(int const k, int const j, int const i) const
    {
        return val;
    }


    L zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};




template<typename L, typename R>
class Expr<L, R, opSum>
{
    const L l;
    const R r;

    typedef decltype(zeroOf(l)*zeroOf(r)) ResType;

public:
    Expr(const L &_l, const R  &_r): l(_l), r(_r) {};

    ResType Eval(int const j, int const i)              const { return l.Eval(j,i)   + r.Eval(j,i); };
    ResType Eval(int const k, int const j, int const i) const { return l.Eval(k,j,i) + r.Eval(k,j,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, typename R>
class Expr<L, R, opSub>
{
    const L l;
    const R r;

    typedef decltype(zeroOf(l)*zeroOf(r)) ResType;

public:
    Expr(const L &_l, const R  &_r): l(_l), r(_r) {};

    ResType Eval(int const j, int const i)              const { return l.Eval(j,i) - r.Eval(j,i); };
    ResType Eval(int const k, int const j, int const i) const { return l.Eval(k,j,i) - r.Eval(k,j,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, typename R>
class Expr<L, R, opMul>
{
    const L l;
    const R r;
    typedef decltype(zeroOf(l)*zeroOf(r)) ResType;

public:
    Expr(const L &_l, const R  &_r): l(_l), r(_r) {};

    ResType Eval(int const j, int const i)              const { return l.Eval(j,i)   * r.Eval(j,i); };
    ResType Eval(int const k, int const j, int const i) const { return l.Eval(k,j,i) * r.Eval(k,j,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, typename R>
class Expr<L, R, opDiv>
{
    const L l;
    const R r;
    typedef decltype(zeroOf(l)*zeroOf(r)) ResType;

public:
    Expr(const L &_l, const R  &_r): l(_l), r(_r) {};

    ResType Eval(int const j, int const i)              const { return l.Eval(j,i)   / r.Eval(j,i); };
    ResType Eval(int const k, int const j, int const i) const { return l.Eval(k,j,i) + r.Eval(k,j,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
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
auto UtoR(const T &U)
{
    return Expr<T, void, opUtoR>(U);
}

template<typename T>
auto VtoR(const T &V)
{
    return Expr<T, void, opVtoR>(V);
}

template<typename T>
auto UtoP(const T &U)
{
    return Expr<T, void, opUtoP>(U);
}

template<typename T>
auto VtoP(const T &V)
{
    return Expr<T, void, opVtoP>(V);
}

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
auto UtoUW_4th(const T &U) {
    return Expr<T, void, opUtoUW_4th>(U);
}


template<typename L, typename R>
auto upwindUtoR(const L &U, const R &velU) {
    return Expr<L, R, opUpwindUtoR>(U, velU);
}

template<typename L, typename R>
auto upwindVtoR(const L &V, const R &velV) {
    return Expr<L, R, opUpwindVtoR>(V, velV);
}

template<typename L, typename R>
auto upwindUtoP(const L &U, const R &velU) {
    return Expr<L, R, opUpwindUtoR>(U, velU);
}

template<typename L, typename R>
auto upwindVtoP(const L &V, const R &velV) {
    return Expr<L, R, opUpwindVtoR>(V, velV);
}




template<typename T>
auto DxRtoU(const T &R, const Expr<double, void, opStencil2D> &on_u)
{
    return Expr<T, void, opDxRtoU>(R, on_u);
}

template<typename T>
auto DyRtoV(const T &R, const Expr<double, void, opStencil2D> &om_v)
{
    return Expr<T, void, opDyRtoV>(R, om_v);
}


template<typename T1, typename T2>
auto divUVtoR(const T1 &U, const T2 &V, const Expr<double, void, opStencil2D> &on_u, const Expr<double, void, opStencil2D> &om_v, const Expr<double, void, opStencil2D> &pn, const Expr<double, void, opStencil2D> &pm)
{
    return Expr<T1, T2, opdivUVtoR>(U, V, on_u, om_v, pn, pm);
}


template<typename T>
auto DXXUtoU(const T &U) { return Expr<T, void, opDXXcentered>(U); }

template<typename T>
auto DXXVtoV(const T &V) { return Expr<T, void, opDXXcentered>(V); }

template<typename T>
auto DEEUtoU(const T &U) { return Expr<T, void, opDEEcentered>(U); }

template<typename T>
auto DEEVtoV(const T &V) { return Expr<T, void, opDEEcentered>(V); }



template<typename T>
auto DXRtoU(const T &R) { return Expr<T, void, opDX>(R); }

template<typename T>
auto DEPtoU(const T &P) { return Expr<T, void, opDE>(P); }


template<typename T>
auto DXPtoV(const T &P) { return Expr<T, void, opDX>(P); }

template<typename T>
auto DERtoV(const T &R) { return Expr<T, void, opDE>(R); }









