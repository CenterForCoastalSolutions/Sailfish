

#define STENCIL(var)      class Expr<double, void, opStencil2D>  const var((double *)(_##var + i))
#define STENCIL3D(var, K) class Expr<double, void, opStencil3D>  const var((double *)(_##var + i), K)

#define STENCILR(var)      class Expr<double, void, opStencil2D, ntR>  const var((double *)(_##var + i))
#define STENCILR3D(var, K) class Expr<double, void, opStencil3D, ntR>  const var((double *)(_##var + i), K)

#define STENCILU(var)      class Expr<double, void, opStencil2D, ntU>  const var((double *)(_##var + i))
#define STENCILU3D(var, K) class Expr<double, void, opStencil3D, ntU>  const var((double *)(_##var + i), K)

#define STENCILV(var)      class Expr<double, void, opStencil2D, ntV>  const var((double *)(_##var + i))
#define STENCILV3D(var, K) class Expr<double, void, opStencil3D, ntV>  const var((double *)(_##var + i), K)

#define STENCILW3D(var, K) class Expr<double, void, opStencil3D, ntW>  const var((double *)(_##var + i), K)





#define OPERATOR_ACCESS ResType operator()(int const k, int const j, int const i) const { return Eval(k,j,i); } \
                        ResType operator()(int const j, int const i) const { return Eval(j,i); }


//#define TO_ACCESS   auto toU() const { if constexpr(nt == ntR) return RtoU(*this); else return 0; } \
//                    auto toV() const { if constexpr(nt == ntR) return RtoV(*this); else return 0; } \
//                    template<NodeType nodeType2> auto to();  // nodeType2 because there is a nodeType in the base template.                                         \
//                    template<> decltype(RtoU(*this)) to<ntU>() const { if constexpr(nt == ntR) return RtoU(*this);  }                               \
//                    template<> decltype(RtoV(*this)) to<ntV>() const { if constexpr(nt == ntR) return RtoV(*this);  }


#define TO_ACCESS   template<NodeType ntTo> auto to() const { if constexpr((nt == ntR) && (ntTo == ntU))  return RtoU(*this); \
                                                              else if constexpr((nt == ntR) && (ntTo == ntV))  return RtoV(*this); \
                                                              else if constexpr((nt == ntW) && (ntTo == ntUW)) return WtoUW(*this); \
                                                              else if constexpr((nt == ntW) && (ntTo == ntVW)) return WtoVW(*this); \
                                                              else {printf("Error Here: %i   %i\n ", int(nt), int(ntTo)); return Expr<double, void, opStencil2D, ntAny>  ((double *)(0));} }
//                                                              static_assert(false, "Member function to() undefined for this choice of nodes");  }


int sz2D, sz3D, szI, szJ, szK;
double *_on_u, *_om_v, *_pm, *_pn;

enum Op  {opValue, opStencil2D, opStencil3D, opUnaryMinus, opSum, opSub, opMul, opDiv, opRtoU, opRtoV, opWtoUW, opWtoVW,
          opUtoR, opVtoR, opUtoP, opVtoP, opUpwindUtoR, opUpwindVtoR, opUpwindUtoP, opUpwindVtoP,
          opRtoU_4th, opRtoV_4th, opRtoW_4th,
          opDxRtoU, opDyRtoV, opdivUVtoR,
          opDXUtoR, opDEVtoR,
          opDXXcentered, opDEEcentered, opDX, opDE,
          opDsigWtoR, opDsigRtoW };

enum NodeType {ntR = 1, ntU = 2, ntV = 4, ntW = 8, ntP = 16, ntUW = 2+8, ntVW = 4+8, ntAny=0, ntInvalid=32};


constexpr bool checkNodeTypes(NodeType const n1, NodeType const n2) { return (n1 == ntAny) || (n2 == ntAny) || (n1 == n2); };
constexpr NodeType compNodes(NodeType const nt1, NodeType const nt2)
{
    if      (nt1 == ntAny) return nt2;
    else if (nt2 == ntAny) return nt1;
    else if (nt1 == nt2)   return nt1;
    else return ntInvalid;
}



//-------------------------------------------


// REMEMBER: This code has to be syncronized with the same definitions in Python code.
constexpr unsigned int bcNone               = 0;
constexpr unsigned int bcAcquire            = 1 << 0;       // Read lateral boundary data from files
constexpr unsigned int bcChapman_explicit   = 1 << 1;
constexpr unsigned int bcChapman_implicit   = 1 << 2;
constexpr unsigned int bcClamped            = 1 << 3;
constexpr unsigned int bcClosed             = 1 << 4;
constexpr unsigned int bcFlather            = 1 << 5;       // (2D momentum)
constexpr unsigned int bcGradient           = 1 << 6;
constexpr unsigned int bcMixed              = 1 << 7;
constexpr unsigned int bcNested             = 1 << 8;
constexpr unsigned int bcNudging            = 1 << 9;
constexpr unsigned int bcPeriodic           = 1 << 10;
constexpr unsigned int bcRadiation          = 1 << 11;
constexpr unsigned int bcReduced            = 1 << 12;      // (2D momentum)
constexpr unsigned int bcShchepetkin        = 1 << 13;      // (2D momentum)


//-------------------------------------------


struct VerticalVelEq
{
    double *D, *SD, *RHS;

    VerticalVelEq(double * _D, double * _SD, double * _RHS): D(_D), SD(_SD), RHS(_RHS) {};
    VerticalVelEq(void) {};

};

//-------------------------------------------


bool isInnerCell(const int idx)
// It is a little unfortunate that i sometimes mean the flat index and other times on the the 2D/3D indices i,j,k.
// In these functions the flat index is called idx to avoid confusion.
// A U node here means that it is a U node not in the boundaries, that is, only an internal node (horizontally).
{
    if (idx >= sz2D) return false;

    const int i = idx % szI;
    const int j = idx / szI;

    return (i > 1) && (i < szI-2 ) && (j > 1) && (j < szJ-2);
}


bool isRNode(const int idx)
// It is a little unfortunate that i sometimes mean the flat index and other times on the the 2D/3D indices i,j,k.
// In these functions the flat index is called idx to avoid confusion.
// An R node here means that it is an R node not in the boundaries, that is, only an internal node (horizontally).
{
    if (idx >= sz2D) return false;

    const int i = idx % szI;
    const int j = idx / szI;

    return (i > 0) && (i < szI) && (j > 0) && (j < szJ);
}

bool isUNode(const int idx)
// It is a little unfortunate that i sometimes mean the flat index and other times one of the 2D/3D indices i,j,k.
// In these functions the flat index is called idx to avoid confusion.
// A U node here means that it is a U node not in the boundaries, that is, only an internal node (horizontally).
{
    if (idx >= sz2D) return false;

    const int i = idx % szI;
    const int j = idx / szI;

    return (i > 1) && (i < szI ) && (j > 0) && (j < szJ);
}

bool isVNode(const int idx)
// It is a little unfortunate that i sometimes mean the flat index and other times on the the 2D/3D indices i,j,k.
// In these functions the flat index is called idx to avoid confusion.
// A V node here means that it is a V node not in the boundaries, that is, only an internal node (horizontally).
{
    if (idx >= sz2D) return false;

    const int i = idx % szI;
    const int j = idx / szI;

    return (i > 0) && (i < szI) && (j > 1) && (j < szJ - 1);
}

extern "C" __global__ void initialize(unsigned int sizeK, unsigned int sizeJ, unsigned int sizeI,
                                      double *on_u, double *om_v, double *pn, double *pm)
{
    printf("Initializing kernel information: %i %i %i \n", sizeK, sizeJ, sizeI);

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



template<class T>
auto zeroOf(T const &t) { return t.zero(); };


template<>
auto zeroOf(double const &t) { return 0.0; };




template<typename L, typename R, Op op, NodeType nodeType = ntAny>
class Expr
{
};


template<typename L, NodeType nodeType>
class Expr<L, void, opStencil2D, nodeType>
{
    L * const p;

public:
    typedef L ResType;
    static const NodeType nt = nodeType;

    Expr(L * _p): p(_p)
    {
    }

    TO_ACCESS

    L &operator()(int const j, int const i) const
    {
        return *(p + j*szI + i);
    }

    // By design, we choose not to return a value from the assignment operator.
    void operator=(L val) const
    {
        *p = val;
    }

    void operator-=(L val) const
    {
        *p -= val;
    }

    void operator+=(L val) const
    {
        *p += val;
    }

    // By design, we choose not to return a value from the assignment operator.
//    template<typename T1>
//    void operator=(const T1 &expr) const
//    {
//        *p = expr.Eval(0,0);
//    }


    template<typename T1>
    void operator=(const T1 &expr) const
    {
        *p = expr.Eval(0,0);
    }


    template<typename T1>
    void operator-=(const T1 &expr) const
    {
        *p = *p - expr.Eval(0,0);
    }

    template<typename T1>
    void operator+=(const T1 &expr) const
    {
        *p = *p + expr.Eval(0,0);
    }


    L Eval(int const j, int const i) const
    {
        return *(p + j*szI + i);
    }

    L Eval(int const k, int const j, int const i) const
    // This function is for compatibility with 3D stencils (for example for operations between 2D and 3D). k is unused.
    {
        return *(p + j*szI + i);
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;

};



template<typename L, NodeType nodeType>
class Expr<L, void, opStencil3D, nodeType>
{
    L * const p;
    int &K;


public:
    static const NodeType nt = nodeType;
    typedef L ResType;

    Expr(L * _p, int &_K): p(_p), K(_K)
    {
    }

    TO_ACCESS

    auto operator[](int const k) const
    {
        Expr<double, void, opStencil2D>  res(p + k*sz2D);
        return res;
    }

    L operator()(int const k, int const j, int const i) const
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

    void operator-=(L val) const
    {
        *(p + K*sz2D) = *(p + K*sz2D) - val;
    }

    template<typename T1>
    void operator-=(const T1 &expr) const
    {
        *(p + K*sz2D) = *(p + K*sz2D) - expr.Eval(0,0,0);
    }

    void operator+=(L val) const
    {
        *(p + K*sz2D) = *(p + K*sz2D) + val;
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
class Expr<L, void, opUtoR, ntR>
{
    L &expr;
    static_assert(L::nt == ntU || L::nt == ntAny);

public:
    static const NodeType nt = ntR;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;


    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS
    TO_ACCESS

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
class Expr<L, void, opVtoR, ntR>
{
    L &expr;
    static_assert(L::nt == ntV || L::nt == ntAny);

public:
    static const NodeType nt = ntR;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;


    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS
    TO_ACCESS


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
class Expr<L, void, opUtoP, ntP>
{
    L &expr;
    static_assert(L::nt == ntU || L::nt == ntAny);

public:
    static const NodeType nt = ntP;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;


    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS
    TO_ACCESS

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
class Expr<L, void, opVtoP, ntP>
{
    L &expr;
    static_assert(L::nt == ntV || L::nt == ntAny);

public:
    static const NodeType nt = ntP;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;


    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS
    TO_ACCESS

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
class Expr<L, void, opRtoU, ntU>
{
    L &expr;
    static_assert((L::nt == ntR) || (L::nt == ntAny));

public:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;
    static const NodeType nt = ntU;

    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS
    TO_ACCESS

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
class Expr<L, void, opWtoUW, ntUW>
{
    L &expr;
    static_assert((L::nt == ntW) || (L::nt == ntAny));

public:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;
    static const NodeType nt = ntUW;

    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS
    TO_ACCESS

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
class Expr<L, void, opWtoVW, ntVW>
{
    L &expr;
    static_assert((L::nt == ntW) || (L::nt == ntAny));

public:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;
    static const NodeType nt = ntVW;

    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS
    TO_ACCESS

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
class Expr<L, void, opRtoV, ntV>
{
    L &expr;
    static_assert((L::nt == ntR) || (L::nt == ntAny));

public:
    static const NodeType nt = ntV;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

    Expr(L _expr): expr(_expr)
    {
    }

    OPERATOR_ACCESS
    TO_ACCESS

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
class Expr<L, R, opUpwindUtoR, ntR>
{
    L &expr;
    R &vel;
    static_assert(L::nt == ntU || L::nt == ntAny);

public:
    static const NodeType nt = ntR;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

    TO_ACCESS

    Expr(L _expr, R _vel): expr(_expr), vel(_vel)
    {
    }

    ResType Eval(int const k, int const j, int const i) const { return (vel.Eval(k,j,i) + vel.Eval(k,j,i+1) >= 0)?expr.Eval(k,j,i):expr.Eval(k,j,i+1); };
    ResType Eval(int const j, int const i) const              { return (vel.Eval(j,i)   + vel.Eval(j,i+1)   >= 0)?expr.Eval(j,i)  :expr.Eval(j,i+1); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, typename R>
class Expr<L, R, opUpwindVtoR, ntR>
{
    L &expr;
    R &vel;
    static_assert(L::nt == ntV || L::nt == ntAny);

public:
    static const NodeType nt = ntR;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

    Expr(L _expr, R _vel): expr(_expr), vel(_vel)
    {
    }

    TO_ACCESS

    ResType Eval(int const k, int const j, int const i) const { return (vel.Eval(k,j,i) + vel.Eval(k,j+1,i) >= 0)?expr.Eval(k,j,i):expr.Eval(k,j+1,i); };
    ResType Eval(int const j, int const i) const              { return (vel.Eval(j,i)   + vel.Eval(j+1,i)   >= 0)?expr.Eval(j,i)  :expr.Eval(j+1,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, typename R>
class Expr<L, R, opUpwindUtoP, ntP>
{
    L &expr;
    R &vel;

    static_assert(L::nt == ntU || L::nt == ntAny);

public:
    static const NodeType nt = ntP;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

    Expr(L _expr, R _vel): expr(_expr), vel(_vel)
    {
    }

    TO_ACCESS

    ResType Eval(int const k, int const j, int const i) const { return (vel.Eval(k,j,i) + vel.Eval(k,j-1,i) >= 0)?expr.Eval(k,j,i):expr.Eval(k,j-1,i); };
    ResType Eval(int const j, int const i) const              { return (vel.Eval(j,i)   + vel.Eval(j-1,i)   >= 0)?expr.Eval(j,i)  :expr.Eval(j-1,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, typename R>
class Expr<L, R, opUpwindVtoP, ntP>
{
    L &expr;
    R &vel;

    static_assert(L::nt == ntV || L::nt == ntAny);

public:
    static const NodeType nt = ntR;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

    Expr(L _expr, R _vel): expr(_expr), vel(_vel)
    {
    }

    TO_ACCESS

    ResType Eval(int const k, int const j, int const i) const { return (vel.Eval(k,j,i) + vel.Eval(k,j,i-1) >= 0)?expr.Eval(k,j,i):expr.Eval(k,j,i-1); };
    ResType Eval(int const j, int const i) const              { return (vel.Eval(j,i)   + vel.Eval(j,i-1)   >= 0)?expr.Eval(j,i)  :expr.Eval(j,i-1); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};



template<typename L, typename R, NodeType nodeType>
class Expr<L, R, opRtoW_4th, nodeType>
{
    L &Var;
    static_assert(L::nt == ntW || L::nt == ntAny);
    static_assert(nodeType == ntW || nodeType == ntUW || nodeType == ntVW || nodeType == ntAny);

public:
    static const NodeType nt = nodeType;
    typedef decltype(zeroOf(Var)*zeroOf(Var)) ResType;

    Expr(L _Var): Var(_Var)
    {
    }

    TO_ACCESS

    ResType Eval(int const k, int const j, int const i) const { return (9.0/16.0)*(Var.Eval(k,j,i) + Var.Eval(k+1,j,i)) - (1.0/16.0)*(Var.Eval(k-1,j,i) + Var.Eval(k+2,j,i)); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, typename R, NodeType nodeType>
class Expr<L, R, opRtoV_4th, nodeType>
{
    L &Var;
    static_assert(L::nt == ntR || L::nt == ntAny);

public:
    static const NodeType nt = nodeType;
    typedef decltype(zeroOf(Var)*zeroOf(Var)) ResType;

    Expr(L _Var): Var(_Var)
    {
    }

    TO_ACCESS

    ResType Eval(int const k, int const j, int const i) const { return (9.0/16.0)*(Var.Eval(k,j,i) + Var.Eval(k,j+1,i)) - (1.0/16.0)*(Var.Eval(k,j-1,i) + Var.Eval(k,j+2,i)); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, typename R, NodeType nodeType>
class Expr<L, R, opRtoU_4th, nodeType>
{
    L &Var;
    static_assert(L::nt == ntU || L::nt == ntAny);

public:
    static const NodeType nt = nodeType;
    typedef decltype(zeroOf(Var)*zeroOf(Var)) ResType;

    Expr(L _Var): Var(_Var)
    {
    }

    TO_ACCESS

    ResType Eval(int const k, int const j, int const i) const { return (9.0/16.0)*(Var.Eval(k,j,i) + Var.Eval(k,j,i+1)) - (1.0/16.0)*(Var.Eval(k,j,i-1) + Var.Eval(k,j,i+2)); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};






template<typename L>
class Expr<L, void, opDxRtoU, ntU>
{
    L &expr;
    static_assert(L::nt == ntR || L::nt == ntAny);

public:
    static const NodeType nt = ntU;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;
    typedef Expr<ResType, void, opStencil2D, ntU> OType;
    OType &on_u;

    Expr(L _expr, OType _on_u): expr(_expr), on_u(_on_u)
    {
    }

    TO_ACCESS

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
class Expr<L, void, opDyRtoV, ntV>
{
    L &expr;
    static_assert(L::nt == ntR || L::nt == ntAny);

public:
    static const NodeType nt = ntV;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;
    typedef Expr<ResType, void, opStencil2D, ntV> OType;
    OType &om_v;

    Expr(L _expr, OType _om_v): expr(_expr), om_v(_om_v)
    {
    }

    TO_ACCESS

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
class Expr<L, void, opDXUtoR, ntR>
{
    L &expr;

    static_assert(L::nt == ntU || L::nt == ntAny);

public:
    static const NodeType nt = ntR;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

    Expr(L _expr): expr(_expr)
    {
    }

    TO_ACCESS

    ResType Eval(int const j, int const i) const
    {
        return expr.Eval(j,i+1) - expr.Eval(j,i);
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return expr.Eval(k,j,i+1) - expr.Eval(k,j,i);
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L>
class Expr<L, void, opDEVtoR, ntR>
{
    L &expr;
    static_assert(L::nt == ntV || L::nt == ntAny);

public:
    static const NodeType nt = ntR;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

    Expr(L _expr): expr(_expr)
    {
    }

    TO_ACCESS

    ResType Eval(int const j, int const i) const
    {
        return expr.Eval(j+1,i) - expr.Eval(j,i);
    }

    ResType Eval(int const k, int const j, int const i) const
    {
        return expr.Eval(k,j+1,i) - expr.Eval(k,j,i);
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};




template<typename L, NodeType nodeType>
class Expr<L, void, opDX, nodeType>
{
    L &expr;

protected:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

public:
    static const NodeType nt = nodeType;

    Expr(L _expr): expr(_expr)
    {
    }

    TO_ACCESS

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


template<typename L, NodeType nodeType>
class Expr<L, void, opDE, nodeType>
{
    L &expr;

protected:
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

public:
    static const NodeType nt = nodeType;

    Expr(L _expr): expr(_expr)
    {
    }

    TO_ACCESS

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



template<typename L, NodeType nodeType>
class Expr<L, void, opDsigWtoR, nodeType>
{
    L &expr;

//    static_assert(L::nt == ntW || L::nt == ntAny);

public:
    static const NodeType nt = nodeType;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

    Expr(L _expr): expr(_expr)
    {
    }

    TO_ACCESS

    ResType Eval(int const k, int const j, int const i) const
    {
        return expr.Eval(k+1,j,i) - expr.Eval(k,j,i);
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};

template<typename L>
class Expr<L, void, opDsigRtoW, ntR>
{
    L &expr;

//    static_assert(L::nt == ntR || L::nt == ntAny);

public:
    static const NodeType nt = ntW;
    typedef decltype(zeroOf(expr)*zeroOf(expr)) ResType;

    Expr(L _expr): expr(_expr)
    {
    }

    TO_ACCESS

    ResType Eval(int const k, int const j, int const i) const
    {
        return expr.Eval(k,j,i) - expr.Eval(k-1,j,i);
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};



template<typename L, typename R>
class Expr<L, R, opdivUVtoR, ntR>
{
    L &U;
    R &V;

    static_assert((L::nt == ntU || L::nt == ntAny) && (R::nt == ntV || R::nt == ntAny));

public:
    static const NodeType nt = ntR;
    typedef decltype(zeroOf(U)*zeroOf(V)) ResType;
    typedef Expr<ResType, void, opStencil2D, ntR> OTypeR;
    typedef Expr<ResType, void, opStencil2D, ntU> OTypeU;
    typedef Expr<ResType, void, opStencil2D, ntV> OTypeV;
    const OTypeU &on_u;
    const OTypeV &om_v;
    const OTypeR &pn;
    const OTypeR &pm;

    Expr(L _U, R _V, OTypeU const &_on_u, OTypeV const &_om_v, OTypeR const &_pn, OTypeR const &_pm): U(_U), V(_V), on_u(_on_u), om_v(_om_v), pn(_pn), pm(_pm)
    {

    }

    TO_ACCESS

    ResType Eval(int const j, int const i) const
    {

        return  ( (U.Eval(j, i+1)*on_u.Eval(j, i+1) - U.Eval(j, i)*on_u.Eval(j, i)) + (V.Eval(j+1, i)*om_v.Eval(j+1, i) - V.Eval(j, i)*om_v.Eval(j, i)) )*pm.Eval(j, i)*pn.Eval(j, i);
    }

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, NodeType nodeType>
class Expr<L, void, opDXXcentered, nodeType>
{
    L &Var;

public:
    static const NodeType nt = L::nt;
    typedef decltype(zeroOf(Var)*zeroOf(Var)) ResType;

    Expr(L _Var): Var(_Var)
    {
    }

    TO_ACCESS

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


template<typename L, NodeType nodeType>
class Expr<L, void, opDEEcentered, nodeType>
{
    L &Var;

public:
    static const NodeType nt = L::nt;
    typedef decltype(zeroOf(Var)*zeroOf(Var)) ResType;

    Expr(L _Var): Var(_Var)
    {
    }

    TO_ACCESS

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
    static const NodeType nt = ntAny;

    Expr(const L _val): val(_val)  {}

    Expr(L &_val): val(_val)  {}

    TO_ACCESS

    const L Eval(int const j, int const i) const
    {
        return val;
    }

    const L Eval(int const k, int const j, int const i) const
    {
        return val;
    }


    L zero(void) const { return L(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};




template<typename L, typename R, NodeType nodeType>
class Expr<L, R, opSum, nodeType>
{
    const L l;
    const R r;

    typedef decltype(zeroOf(l)*zeroOf(r)) ResType;
    static_assert(checkNodeTypes(L::nt,R::nt));

public:
    static const NodeType nt = nodeType;

    Expr(const L &_l, const R  &_r): l(_l), r(_r) {};

    TO_ACCESS

    ResType Eval(int const j, int const i)              const { return l.Eval(j,i)   + r.Eval(j,i); };
    ResType Eval(int const k, int const j, int const i) const { return l.Eval(k,j,i) + r.Eval(k,j,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};




template<typename L, NodeType nodeType>
class Expr<L, void, opUnaryMinus, nodeType>
{
    const L l;

    typedef decltype(zeroOf(l)) ResType;

public:
    static const NodeType nt = nodeType;

    Expr(const L &_l): l(_l) {};

    TO_ACCESS

    ResType Eval(int const j, int const i)              const { return -l.Eval(j,i); };
    ResType Eval(int const k, int const j, int const i) const { return -l.Eval(k,j,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, typename R, NodeType nodeType>
class Expr<L, R, opSub, nodeType>
{
    const L l;
    const R r;

    typedef decltype(zeroOf(l)*zeroOf(r)) ResType;
    static_assert(checkNodeTypes(L::nt, R::nt));

public:
    static const NodeType nt = nodeType;

    Expr(const L &_l, const R  &_r): l(_l), r(_r) {};

    TO_ACCESS

    ResType Eval(int const j, int const i)              const { return l.Eval(j,i)   - r.Eval(j,i); };
    ResType Eval(int const k, int const j, int const i) const { return l.Eval(k,j,i) - r.Eval(k,j,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, typename R, NodeType nodeType>
class Expr<L, R, opMul, nodeType>
{
    const L l;
    const R r;
    typedef decltype(zeroOf(l)*zeroOf(r)) ResType;
    static_assert(checkNodeTypes(L::nt, R::nt));

public:
    static const NodeType nt = nodeType;

    Expr(const L &_l, const R  &_r): l(_l), r(_r) {};

    TO_ACCESS

    ResType Eval(int const j, int const i)              const { return l.Eval(j,i)   * r.Eval(j,i); };
    ResType Eval(int const k, int const j, int const i) const { return l.Eval(k,j,i) * r.Eval(k,j,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};


template<typename L, typename R, NodeType nodeType>
class Expr<L, R, opDiv, nodeType>
{
    const L l;
    const R r;
    typedef decltype(zeroOf(l)*zeroOf(r)) ResType;
    static_assert(checkNodeTypes(L::nt, R::nt));

public:
    static const NodeType nt = nodeType;

    Expr(const L &_l, const R  &_r): l(_l), r(_r) {};

    TO_ACCESS

    ResType Eval(int const j, int const i)              const { return l.Eval(j,i)   / r.Eval(j,i); };
    ResType Eval(int const k, int const j, int const i) const { return l.Eval(k,j,i) / r.Eval(k,j,i); };

    ResType zero(void) const { return ResType(0.0); }  // This is mostly a trick to get a number with the appropriate type;
};







template<typename L>
auto operator-(const L &l) { return Expr<L, void, opUnaryMinus>(l); }

template<typename L, typename R>
auto operator-(const L &l, const R &r) { return Expr<L, R, opSub, compNodes(L::nt, R::nt)>(l, r); }

template<typename L, typename R>
auto operator*(const L &l, const R &r) { return Expr<L, R, opMul, compNodes(L::nt, R::nt)>(l, r); }

template<typename L, typename R>
auto operator/(const L &l, const R &r) { return Expr<L, R, opDiv, compNodes(L::nt, R::nt)>(l, r); }

template<typename L, typename R>
auto operator+(const L &l, const R &r) { return Expr<L, R, opSum, compNodes(L::nt, R::nt)>(l, r); }




template<typename R>
auto operator+(const double l, const R &r) { auto ll = Expr<double, void, opValue>(l); return Expr<decltype(ll), R, opSum, R::nt>(ll, r); }

template<typename R>
auto operator-(const double l, const R &r) { auto ll = Expr<double, void, opValue>(l); return Expr<decltype(ll), R, opSub, R::nt>(ll, r); }

template<typename R>
auto operator*(const double l, const R &r) { auto ll = Expr<double, void, opValue>(l); return Expr<decltype(ll), R, opMul, R::nt>(ll, r); }

template<typename R>
auto operator/(const double l, const R &r) { auto ll = Expr<double, void, opValue>(l); return Expr<decltype(ll), R, opDiv, R::nt>(ll, r); }




template<typename L>
auto operator+(const L &l, const double r) { auto rr = Expr<double, void, opValue>(r); return Expr<L, decltype(rr), opSum, L::nt>(l, rr); }

template<typename L>
auto operator-(const L &l, const double r) { auto rr = Expr<double, void, opValue>(r); return Expr<L, decltype(rr), opSub, L::nt>(l, rr); }

template<typename L>
auto operator*(const L &l, const double r) { auto rr = Expr<double, void, opValue>(r); return Expr<L, decltype(rr), opMul, L::nt>(l, rr); }

template<typename L>
auto operator/(const L &l, const double r) { auto rr = Expr<double, void, opValue>(r); return Expr<L, decltype(rr), opDiv, L::nt>(l, rr); }


template<typename T>
auto UtoR(const T &U)
{
    return Expr<T, void, opUtoR, ntR>(U);
}

template<typename T>
auto VtoR(const T &V)
{
    return Expr<T, void, opVtoR, ntR>(V);
}

template<typename T>
auto UtoP(const T &U)
{
    return Expr<T, void, opUtoP, ntP>(U);
}

template<typename T>
auto VtoP(const T &V)
{
    return Expr<T, void, opVtoP, ntP>(V);
}

template<typename T>
auto RtoU(const T &R)
{
    return Expr<T, void, opRtoU, ntU>(R);
}


template<typename T>
auto RtoV(const T &R)
{
    return Expr<T, void, opRtoV, ntV>(R);
}


template<typename T>
auto WtoUW(const T &W) {
    return Expr<T, void, opWtoUW, ntUW>(W);
}

template<typename T>
auto WtoVW(const T &W) {
    return Expr<T, void, opWtoVW, ntVW>(W);
}



template<typename T>
auto RtoW_4th(const T &R) {
    return Expr<T, void, opRtoW_4th, ntW>(R);
}

template<typename T>
auto UtoUW_4th(const T &U) {
    return Expr<T, void, opRtoW_4th, ntUW>(U);
}

template<typename T>
auto VtoVW_4th(const T &V) {
    return Expr<T, void, opRtoW_4th, ntVW>(V);
}

template<typename T>
auto WtoUW_4th(const T &W) {
    return Expr<T, void, opRtoU_4th, ntUW>(W);
}

template<typename T>
auto WtoVW_4th(const T &W) {
    return Expr<T, void, opRtoV_4th, ntVW>(W);
}


template<typename L, typename R>
auto upwindUtoR(const L &U, const R &velU) {
    return Expr<L, R, opUpwindUtoR, ntR>(U, velU);
}

template<typename L, typename R>
auto upwindVtoR(const L &V, const R &velV) {
    return Expr<L, R, opUpwindVtoR, ntR>(V, velV);
}

template<typename L, typename R>
auto upwindUtoP(const L &U, const R &velU) {
    return Expr<L, R, opUpwindUtoP, ntP>(U, velU);  // TODO: CHECK if this and the next are correct (opUpwindUtoR instead of opUpwindUtoP).
}

template<typename L, typename R>
auto upwindVtoP(const L &V, const R &velV) {
    return Expr<L, R, opUpwindVtoP, ntP>(V, velV);
}




template<typename T>
auto DXUtoR(const T &R)
{
    return Expr<T, void, opDXUtoR, ntR>(R);
}

template<typename T>
auto DEVtoR(const T &R)
{
    return Expr<T, void, opDEVtoR, ntR>(R);
}

template<typename T>
auto DxRtoU(const T &R, const Expr<double, void, opStencil2D> &on_u)
{
    return Expr<T, void, opDxRtoU, ntU>(R, on_u);
}

template<typename T>
auto DyRtoV(const T &R, const Expr<double, void, opStencil2D> &om_v)
{
    return Expr<T, void, opDyRtoV, ntV>(R, om_v);
}




typedef class Expr<double, void, opStencil2D, ntR> DblStencil2DR;
typedef class Expr<double, void, opStencil2D, ntU> DblStencil2DU;
typedef class Expr<double, void, opStencil2D, ntV> DblStencil2DV;
//, const Expr<double, void, opStencil2D, ntU> &on_u, const Expr<double, void, opStencil2D, ntV> &om_v, const Expr<double, void, opStencil2D, ntR> &pn, const Expr<double, void, opStencil2D, ntR> &pm    , on_u, om_v, pn, pm
template<typename T1, typename T2>
auto divUVtoR(const T1 &U, const T2 &V, DblStencil2DU const &_on_u, DblStencil2DV const &_om_v, DblStencil2DR const &_pn, DblStencil2DR const &_pm)
{
    return Expr<T1, T2, opdivUVtoR, ntR>(U, V, _on_u, _om_v, _pn, _pm);
}





template<typename T>
auto DXXUtoU(const T &U) { return Expr<T, void, opDXXcentered, ntU>(U); }

template<typename T>
auto DXXVtoV(const T &V) { return Expr<T, void, opDXXcentered, ntV>(V); }

template<typename T>
auto DEEUtoU(const T &U) { return Expr<T, void, opDEEcentered, ntU>(U); }

template<typename T>
auto DEEVtoV(const T &V) { return Expr<T, void, opDEEcentered, ntV>(V); }



template<typename T>
auto DXRtoU(const T &R) { return Expr<T, void, opDX, ntU>(R); }

template<typename T>
auto DEPtoU(const T &P) { return Expr<T, void, opDE, ntU>(P); }


template<typename T>
auto DXPtoV(const T &P) { return Expr<T, void, opDX, ntV>(P); }

template<typename T>
auto DERtoV(const T &R) { return Expr<T, void, opDE, ntV>(R); }



template<typename T>
auto DsigUWtoU(const T &UW) { return Expr<T, void, opDsigWtoR, ntUW>(UW); }


template<typename T>
auto DsigVWtoV(const T &UW) { return Expr<T, void, opDsigWtoR, ntV>(UW); }


template<typename T>
auto DsigWtoR(const T &W) { return Expr<T, void, opDsigWtoR, ntR>(W); }


template<typename T>
auto DsigRtoW(const T &W) { return Expr<T, void, opDsigRtoW, ntR>(W); }


// ---------------------------------------------------------------


