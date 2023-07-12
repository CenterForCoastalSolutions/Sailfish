enum Op  {opValue, opStencil, opSum, opSub, opMul, opDiv, opRtoU, opRtoV,
          opDERtoU, opDNRtoV, opdivUVtoR};


extern "C" __global__ void initialize(unsigned int sizeK, unsigned int sizeJ, unsigned int sizeI,
                                      double *_on_u, double *_om_v, double *_pm, double *_pn);

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

    // By design, we choose not to return a value from the assignment opertor.
    void operator=(L val) const
    {
        *p = val;
    }

    // By design, we choose not to return a value from the assignment opertor.
    template<typename T1>
    void operator=(const T1 &expr) const
    {
        *p = expr.Eval(0,0);
    }


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

public:
    typedef decltype(expr.Eval(0,0)*expr.Eval(0,0)) ResType;

    Expr(L _expr): expr(_expr)
    {
    }

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

