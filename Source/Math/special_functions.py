# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "marimo>=0.19.6",
#     "numpy==2.4.1",
#     "scipy==1.17.0",
#     "sympy==1.14.0",
# ]
# ///

import marimo

__generated_with = "0.19.6"
app = marimo.App(width="medium")


@app.cell
def _():
    import numpy as np
    import sympy as sp
    from sympy.codegen.ast import real, float32, float64
    from sympy.codegen.rewriting import optimize, optims_c99
    import marimo as mo
    N_MAX = 16
    case_indent = '\t\t'
    code_indent = '\t\t\t'
    return (
        N_MAX,
        case_indent,
        code_indent,
        float32,
        float64,
        optimize,
        optims_c99,
        real,
        sp,
    )


@app.cell
def _(case_indent, code_indent, float32, optimize, optims_c99, real, sp):
    def generate_Laguerre(n: int, device_library='', type=float32):
        if type == float32:
            zero = '0.0f'
        else:
            zero = '0.0'

        x = sp.symbols('x')
        a = sp.symbols('alpha')
        core = sp.exp(-x) * x**(n+a)
        prefactor = x**(-a) * sp.exp(x) / sp.factorial(n)
        rodrigues =  optimize(sp.simplify(prefactor * core.diff(x, n)), optimizations=optims_c99)
        code = sp.cxxcode(rodrigues, type_aliases={real: type})
        statement = f'return {(code)};'
        return code_indent + statement.replace('std::',device_library).replace('f(','(') + '\n'

    def generate_Hermite(n: int, device_library='', type=float32):
        x = sp.symbols('x')
        rodrigues =  optimize(sp.simplify((((-1)**n) * sp.exp(x**2)) * (sp.exp(-x**2).diff(x, n))), optimizations=optims_c99)
        code = sp.cxxcode(rodrigues, type_aliases={real: type})
        statement = f'return {code};'
        return code_indent + statement.replace('std::',device_library).replace('f(','(') + '\n'

    def generate_case(n: int, case_code: str):
        return case_indent + f'case {n}u:' + '\n' + case_code# + '\n' + code_indent + 'break;' '\n'
    return generate_Hermite, generate_Laguerre, generate_case


@app.cell
def _(
    N_MAX,
    case_indent,
    float64,
    generate_Hermite,
    generate_Laguerre,
    generate_case,
):
    laguerre_code = r"""deviceFunction static inline constexpr float generalized_laguerre(const float x, const int alpha, const uint8_t n){
        switch(n){
    """
    for i in range(N_MAX):
        laguerre_code += generate_case(i,generate_Laguerre(i,device_library='deviceFPLib::'))
    laguerre_code += case_indent + 'default: return 0.0f;\n'
    laguerre_code += r"""    }
    }
    """

    laguerre_code_double = r"""deviceFunction static inline constexpr double generalized_laguerre(const double x, const int alpha, const uint8_t n){
        switch(n){
    """
    for i in range(N_MAX):
        laguerre_code_double += generate_case(i,generate_Laguerre(i,device_library='deviceFPLib::',type=float64))
    laguerre_code_double += case_indent + 'default: return 0.0;\n'
    laguerre_code_double += r"""    }
    }
    """

    hermite_code = r"""deviceFunction static inline constexpr float hermite(const float x, const uint8_t n){
        switch(n){
    """
    for i in range(N_MAX):
        hermite_code += generate_case(i,generate_Hermite(i,device_library='deviceFPLib::'))
    hermite_code += case_indent + 'default: return 0.0f;\n'
    hermite_code += r"""    }
    }
    """

    hermite_code_double = r"""deviceFunction static inline constexpr double hermite(const double x, const uint8_t n){
        switch(n){
    """
    for i in range(N_MAX):
        hermite_code_double += generate_case(i,generate_Hermite(i,device_library='deviceFPLib::', type=float64))
    hermite_code_double += case_indent + 'default: return 0.0;\n'
    hermite_code_double += r"""    }
    }
    """
    return (
        hermite_code,
        hermite_code_double,
        laguerre_code,
        laguerre_code_double,
    )


@app.cell
def _(hermite_code, hermite_code_double, laguerre_code, laguerre_code_double):
    hermite = hermite_code + hermite_code_double
    laguerre = laguerre_code + laguerre_code_double
    return hermite, laguerre


@app.cell
def _(hermite, laguerre):
    with open('hermite.hpp', 'w') as f:
        f.write(hermite)

    with open('laguerre.hpp', 'w') as f:
        f.write(laguerre)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
