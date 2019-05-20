from sage.all import show, assume, diff, sqrt, exp, latex, log, prod, vector, matrix, identity_matrix
from sage.all import var as SageVariable
from sage.all import SR as Symbols, RR as Reals, ZZ as Integers
from sage.all import SageObject

from sys import exit

from poblano import PoblanoRing, PoblanoExpression


def test_constant_poblano_expression_diff_and_jump():
    s = PoblanoRing(['x', 'y'])
    x, y = s.vars
    
    a = s.fxn('a')
    b = s.fxn('b')
    
    return bool(s.jump_coeff(a, x) == 0) and bool(s.jump_coeff(a, y) == 0) and bool(s.jump_coeff(b, x) == 0) and bool(s.jump_coeff(b, y) == 0)

def test_poblano_expression_construction_options():
    s = PoblanoRing(['x', 'y'])
    x, y = s.vars
    d = s.d
    
    jfx = var('jfx')
    jfy = var('jfy')
    
    a = s.fxn('a')
    b = s.fxn('b')
    fv = s.fxn('fv', value=a * x + b * y)
    fd = s.fxn('fd', total_diff=a*d(x) + b*d(y))
    fj = s.fxn('fj', jump_dict={x: a, y: b})
    fvd = s.fxn('fvd', value=a * x + b * y, total_diff=a*d(x) + b*d(y))
    fvj = s.fxn('fvj', value=a * x + b * y, jump_dict={x: jfx, y: jfy})
    fdj = s.fxn('fdj', total_diff=a*d(x) + b*d(y), jump_dict={x: a, y: b})
    fvdj = s.fxn('fvdj', value=a * x + b * y, total_diff=a*d(x) + b*d(y), jump_dict={x: a, y: b})
    
    exact_total_diff = fv.value.diff(x) * d(x) + fv.value.diff(y) * d(y)
    exact_jump_expr_dict = {x: a, y: b}
    exact_jump_symb_dict = {x: a.symbol, y: b.symbol}
    zero_jump_dict = {x: 0, y: 0}
    spec_jump_dict = {x: jfx, y: jfy}
    
    test_passed = True
    
    test_passed = test_passed and bool(fv.total_diff == exact_total_diff)
    test_passed = test_passed and bool(fd.total_diff == exact_total_diff)
    test_passed = test_passed and bool(fj.total_diff == exact_total_diff)
    test_passed = test_passed and bool(fvd.total_diff == exact_total_diff)
    test_passed = test_passed and bool(fvj.total_diff == exact_total_diff)
    test_passed = test_passed and bool(fdj.total_diff == exact_total_diff)
    test_passed = test_passed and bool(fvdj.total_diff == exact_total_diff)
    
    test_passed = test_passed and bool(fv.jump_dict == exact_jump_symb_dict)
    test_passed = test_passed and bool(fd.jump_dict == zero_jump_dict)
    test_passed = test_passed and bool(fj.jump_dict == exact_jump_expr_dict)
    test_passed = test_passed and bool(fvd.jump_dict == exact_jump_symb_dict)
    test_passed = test_passed and bool(fvj.jump_dict == spec_jump_dict)
    test_passed = test_passed and bool(fdj.jump_dict == exact_jump_expr_dict)
    test_passed = test_passed and bool(fvdj.jump_dict == exact_jump_expr_dict)

    return test_passed

def test_poblano_expression_vs_space_diff_and_jump():
    s = PoblanoRing(['x', 'y'])
    x, y = s.vars
    d = s.d
    
    jfx = var('jfx')
    jfy = var('jfy')
    
    a = s.fxn('a')
    b = s.fxn('b')
    fv = s.fxn('fv', value=a * x + b * y)
    fd = s.fxn('fd', total_diff=a*d(x) + b*d(y))
    fj = s.fxn('fj', jump_dict={x: a, y: b})
    fvd = s.fxn('fvd', value=a * x + b * y, total_diff=a*d(x) + b*d(y))
    fvj = s.fxn('fvj', value=a * x + b * y, jump_dict={x: jfx, y: jfy})
    fdj = s.fxn('fdj', total_diff=a*d(x) + b*d(y), jump_dict={x: a, y: b})
    fvdj = s.fxn('fvdj', value=a * x + b * y, total_diff=a*d(x) + b*d(y), jump_dict={x: a, y: b})
    
    test_passed = True
    
    test_passed = test_passed and bool(fv.diff(x) == s.diff(fv, x)) and bool(fv.diff(y) == s.diff(fv, y))
    test_passed = test_passed and bool(fd.diff(x) == s.diff(fd, x)) and bool(fd.diff(y) == s.diff(fd, y))
    test_passed = test_passed and bool(fj.diff(x) == s.diff(fj, x)) and bool(fj.diff(y) == s.diff(fj, y))
    test_passed = test_passed and bool(fvd.diff(x) == s.diff(fvd, x)) and bool(fvd.diff(y) == s.diff(fvd, y))
    test_passed = test_passed and bool(fvj.diff(x) == s.diff(fvj, x)) and bool(fvj.diff(y) == s.diff(fvj, y))
    test_passed = test_passed and bool(fdj.diff(x) == s.diff(fdj, x)) and bool(fdj.diff(y) == s.diff(fdj, y))
    test_passed = test_passed and bool(fvdj.diff(x) == s.diff(fvdj, x)) and bool(fvdj.diff(y) == s.diff(fvdj, y))
    
    test_passed = test_passed and bool(fv.jump_coeff(x) == s.jump_coeff(fv, x)) and bool(fv.jump_coeff(y) == s.jump_coeff(fv, y))
    test_passed = test_passed and bool(fd.jump_coeff(x) == s.jump_coeff(fd, x)) and bool(fd.jump_coeff(y) == s.jump_coeff(fd, y))
    test_passed = test_passed and bool(fj.jump_coeff(x) == s.jump_coeff(fj, x)) and bool(fj.jump_coeff(y) == s.jump_coeff(fj, y))
    test_passed = test_passed and bool(fvd.jump_coeff(x) == s.jump_coeff(fvd, x)) and bool(fvd.jump_coeff(y) == s.jump_coeff(fvd, y))
    test_passed = test_passed and bool(fvj.jump_coeff(x) == s.jump_coeff(fvj, x)) and bool(fvj.jump_coeff(y) == s.jump_coeff(fvj, y))
    test_passed = test_passed and bool(fdj.jump_coeff(x) == s.jump_coeff(fdj, x)) and bool(fdj.jump_coeff(y) == s.jump_coeff(fdj, y))
    test_passed = test_passed and bool(fvdj.jump_coeff(x) == s.jump_coeff(fvdj, x)) and bool(fvdj.jump_coeff(y) == s.jump_coeff(fvdj, y))

    return test_passed

def test_jump_expansion_accuracy_simple_poblano_expression():
    s = PoblanoRing(['x', 'y'])
    x, y = s.vars
    d = s.d
    
    h_jx = SageVariable('h_jx', latex_name='\\mathcal{J}^h_x')
    h_jy = SageVariable('h_jy', latex_name='\\mathcal{J}^h_y')
    h = s.fxn('h', jump_dict={x: h_jx, y: h_jy})
    
    return bool(h.jump_coeff(x) == h_jx) and bool(h.jump_coeff(y) == h_jy) and bool(s.jump_coeff(h, x) == h_jx) and bool(s.jump_coeff(h, y) == h_jy)
    
def test_jump_expansion_accuracy_poblano_expression_and_constants():
    s = PoblanoRing(['x', 'y'])
    x, y = s.vars
    d = s.d
    
    h_jx = SageVariable('h_jx', latex_name='\\mathcal{J}^h_x')
    h_jy = SageVariable('h_jy', latex_name='\\mathcal{J}^h_y')
    h = s.fxn('h', jump_dict={x: h_jx, y: h_jy})
    a = s.fxn('a')
    
    test_passed = True
    
    test_passed = test_passed and bool(s.jump_coeff(2 + h, x) == h_jx)
    test_passed = test_passed and bool(s.jump_coeff(2 * h, x) == 2 * h_jx)
    test_passed = test_passed and bool(s.jump_coeff(a + h, x) == h_jx)
    test_passed = test_passed and bool(s.jump_coeff(a * h, x) == a * h_jx)
    test_passed = test_passed and bool(s.jump_coeff(sqrt(2) + h, x) == h_jx)
    test_passed = test_passed and bool(s.jump_coeff(sqrt(2) * h, x) == sqrt(2) * h_jx)
    
    test_passed = test_passed and bool(s.jump_coeff(2 + h, y) == h_jy)
    test_passed = test_passed and bool(s.jump_coeff(2 * h, y) == 2 * h_jy)
    test_passed = test_passed and bool(s.jump_coeff(a + h, y) == h_jy)
    test_passed = test_passed and bool(s.jump_coeff(a * h, y) == a * h_jy)
    test_passed = test_passed and bool(s.jump_coeff(sqrt(2) + h, y) == h_jy)
    test_passed = test_passed and bool(s.jump_coeff(sqrt(2) * h, y) == sqrt(2) * h_jy)
    
    return test_passed
    
def test_jump_expansion_accuracy_binary_operations():
    s = PoblanoRing(['x', 'y'])
    x, y = s.vars
    d = s.d
    
    h_jx = SageVariable('h_jx', latex_name='\\mathcal{J}^h_x')
    h_jy = SageVariable('h_jy', latex_name='\\mathcal{J}^h_y')
    i_jx = SageVariable('i_jx', latex_name='\\mathcal{J}^i_x')
    i_jy = SageVariable('i_jy', latex_name='\\mathcal{J}^i_y')
    j_jx = SageVariable('j_jx', latex_name='\\mathcal{J}^j_x')
    j_jy = SageVariable('j_jy', latex_name='\\mathcal{J}^j_y')
    
    h = s.fxn('h', jump_dict={x: h_jx, y: h_jy})
    i = s.fxn('i', jump_dict={x: i_jx, y: i_jy})
    j = s.fxn('j', jump_dict={x: j_jx, y: j_jy})
    a = s.fxn('a')
    
    test_passed = True
    
    test_passed = test_passed and bool(s.jump_coeff(h + i, x) == h_jx + i_jx)
    test_passed = test_passed and bool(s.jump_coeff(h + i + j, x) == h_jx + i_jx + j_jx)
    test_passed = test_passed and bool(s.jump_coeff(h - i, x) == h_jx - i_jx)
    test_passed = test_passed and bool(s.jump_coeff(h - i - j, x) == h_jx - i_jx - j_jx)
    test_passed = test_passed and bool(s.jump_coeff(-h - i, x) == -h_jx - i_jx)
    test_passed = test_passed and bool(s.jump_coeff(-h - i - j, x) == -h_jx - i_jx - j_jx)
    test_passed = test_passed and bool(s.jump_coeff(h * i, x) == s.avg(i) * h_jx + s.avg(h) * i_jx)
    test_passed = test_passed and bool(s.jump_coeff(h * i * j, x) == s.avg(i * j) * h_jx + s.avg(h) * (s.avg(j) * i_jx + s.avg(i) * j_jx))
    test_passed = test_passed and bool(s.jump_coeff(h / i, x) == s.avg(1 / i) * h_jx + s.avg(h) / s.avg(i, 'inverse') * i_jx)
    test_passed = test_passed and bool(s.jump_coeff(i / h, x) == s.avg(1 / h) * i_jx + s.avg(i) / s.avg(h, 'inverse') * h_jx)
    
    return test_passed
    
def test_jump_expansion_accuracy_integer_powers():
    s = PoblanoRing(['x', 'y'])
    x, y = s.vars
    d = s.d
    
    h_jx = SageVariable('h_jx', latex_name='\\mathcal{J}^h_x')
    h_jy = SageVariable('h_jy', latex_name='\\mathcal{J}^h_y')
    i_jx = SageVariable('i_jx', latex_name='\\mathcal{J}^i_x')
    i_jy = SageVariable('i_jy', latex_name='\\mathcal{J}^i_y')
    j_jx = SageVariable('j_jx', latex_name='\\mathcal{J}^j_x')
    j_jy = SageVariable('j_jy', latex_name='\\mathcal{J}^j_y')
    
    h = s.fxn('h', jump_dict={x: h_jx, y: h_jy})
    i = s.fxn('i', jump_dict={x: i_jx, y: i_jy})
    j = s.fxn('j', jump_dict={x: j_jx, y: j_jy})
    a = s.fxn('a')
    
    test_passed = True
    
    test_passed = test_passed and bool(s.jump_coeff(i * i, x) == 2 * s.avg(i) * i_jx)
    test_passed = test_passed and bool(s.jump_coeff(i ** 2, x) == 2 * s.avg(i) * i_jx)
    test_passed = test_passed and bool(s.jump_coeff(i ** 4, x) == 4 * s.avg(i) * s.avg(i * i) * i_jx)
    test_passed = test_passed and bool(s.jump_coeff(i ** 8, x) == 8 * s.avg(i) * s.avg(i * i) * s.avg(i * i * i * i) * i_jx)
    test_passed = test_passed and bool(s.jump_coeff(i ** 3, x) == s.avg(i) * s.jump_coeff(i ** 2, x) + s.avg(i ** 2) * i_jx)
    test_passed = test_passed and bool(s.jump_coeff(i ** 5, x) == s.avg(i) * s.jump_coeff(i ** 4, x) + s.avg(i ** 4) * i_jx)
    test_passed = test_passed and bool(s.jump_coeff(i ** 6, x) == 2 * s.avg(i ** 3) * s.jump_coeff(i ** 3, x))
    test_passed = test_passed and bool(s.jump_coeff(i ** 7, x) == s.avg(i) * s.jump_coeff(i ** 6, x) + s.avg(i ** 6) * i_jx)
    
    return test_passed
    
def test_jump_expansion_accuracy_noninteger_powers():
    s = PoblanoRing(['x', 'y'])
    x, y = s.vars
    d = s.d
    
    h_jx = SageVariable('h_jx', latex_name='\\mathcal{J}^h_x')
    h_jy = SageVariable('h_jy', latex_name='\\mathcal{J}^h_y')
    i_jx = SageVariable('i_jx', latex_name='\\mathcal{J}^i_x')
    i_jy = SageVariable('i_jy', latex_name='\\mathcal{J}^i_y')
    j_jx = SageVariable('j_jx', latex_name='\\mathcal{J}^j_x')
    j_jy = SageVariable('j_jy', latex_name='\\mathcal{J}^j_y')
    
    h = s.fxn('h', jump_dict={x: h_jx, y: h_jy})
    i = s.fxn('i', jump_dict={x: i_jx, y: i_jy})
    j = s.fxn('j', jump_dict={x: j_jx, y: j_jy})
    a = s.fxn('a')
    
    test_passed = True
    
    test_passed = test_passed and bool(s.jump_coeff(i ** 0.5, x) == i_jx / (2 * s.avg(i ** 0.5)))
    test_passed = test_passed and bool(s.jump_coeff(i ** -0.5, x) == i_jx / (2 * s.avg(i ** 0.5) * s.avg(i ** 0.5, 'inverse')))
    test_passed = test_passed and bool(s.jump_coeff(i ** (1 / 3), x) == i_jx * 1 / 3 / (s.avg(1 / 3 * log(i.symbol), 'exp') * s.avg(i, 'logarithmic')))
    test_passed = test_passed and bool(s.jump_coeff(i ** (-1 / 3), x) == i_jx * -1 / 3 / (s.avg(-1 / 3 * log(i.symbol), 'exp') * s.avg(i, 'logarithmic')))
    
    return test_passed
    
def test_jump_expansion_accuracy_unary_operations():
    s = PoblanoRing(['x', 'y'])
    x, y = s.vars
    d = s.d
    
    h_jx = SageVariable('h_jx', latex_name='\\mathcal{J}^h_x')
    h_jy = SageVariable('h_jy', latex_name='\\mathcal{J}^h_y')
    i_jx = SageVariable('i_jx', latex_name='\\mathcal{J}^i_x')
    i_jy = SageVariable('i_jy', latex_name='\\mathcal{J}^i_y')
    j_jx = SageVariable('j_jx', latex_name='\\mathcal{J}^j_x')
    j_jy = SageVariable('j_jy', latex_name='\\mathcal{J}^j_y')
    
    h = s.fxn('h', jump_dict={x: h_jx, y: h_jy})
    i = s.fxn('i', jump_dict={x: i_jx, y: i_jy})
    j = s.fxn('j', jump_dict={x: j_jx, y: j_jy})
    a = s.fxn('a')
    
    test_passed = True
    
    test_passed = test_passed and bool(s.jump_coeff(log(h.symbol), x) == h_jx / s.avg(h, 'logarithmic'))
    test_passed = test_passed and bool(s.jump_coeff(exp(h.symbol), x) == h_jx / s.avg(h, 'exp'))
    
    return test_passed
    
def test_jump_expansion_accuracy_polynomial_consistency():
    s = PoblanoRing(['x', 'y'])
    x, y = s.vars
    d = s.d
    
    a_m4 = SageVariable('a_m4')
    a_m3 = SageVariable('a_m3')
    a_m2 = SageVariable('a_m2')
    a_m1 = SageVariable('a_m1')
    a_0 = SageVariable('a_0')
    a_1 = SageVariable('a_1')
    a_2 = SageVariable('a_2')
    a_3 = SageVariable('a_3')
    a_4 = SageVariable('a_4')
    a_5 = SageVariable('a_5')
    a_6 = SageVariable('a_6')
    a_7 = SageVariable('a_7')
    poly = a_m4 * x ** -4 + a_m3 * x ** -3 + a_m2 * x ** -2 + a_m1 * x ** -1 + a_0 + a_1 * x + a_2 * x ** 2 + a_3 * x ** 3 + a_4 * x ** 4 + a_5 * x ** 5 + a_6 * x ** 6 + a_7 * x ** 7
    jpx = s.jump_coeff(poly, x)
    return bool(poly.diff(x) == s.subs_for_consistency(jpx))


def test_differentiation_over_constants_and_expands():
    s = PoblanoRing(['x', 'y'], {'x': '\\mathrm{x}'})
    x, y = s.vars
    d = s.d

    a = SageVariable('a')
    b = SageVariable('b')

    af = s.fxn('a_f')
    bf = s.fxn('b_f')

    c = s.fxn('c', 2 * sqrt(2) * af + bf)

    # we test across constants that are numbers (2, 3), numeric symbols (sqrt(2), sqrt(3)), Sage variables (a, b), and constant PoblanoExpression objects (af, bf)
    f_1 = s.fxn('f_1', 2 * x + 3 * y)
    f_2 = s.fxn('f_2', sqrt(2) * x + sqrt(3) * y)
    f_3 = s.fxn('f_3', a * x + b * y)
    f_4 = s.fxn('f_4', af * x + bf * y)

    # g_is is given in a condensed form
    g_1s = s.fxn('g_1s', (x * y * (2 * sqrt(2) * c * f_1 + f_1 ** 2 / y) + exp(f_1 * x)))
    g_2s = s.fxn('g_2s', (x * y * (2 * sqrt(2) * c * f_2 + f_2 ** 2 / y) + exp(f_2 * x)))
    g_3s = s.fxn('g_3s', (x * y * (2 * sqrt(2) * c * f_3 + f_3 ** 2 / y) + exp(f_3 * x)))
    g_4s = s.fxn('g_4s', (x * y * (2 * sqrt(2) * c * f_4 + f_4 ** 2 / y) + exp(f_4 * x)))

    # g_ie is expanded
    g_1e = s.fxn('g_1e', g_1s.value.expand())
    g_2e = s.fxn('g_2e', g_2s.value.expand())
    g_3e = s.fxn('g_3e', g_3s.value.expand())
    g_4e = s.fxn('g_4e', g_4s.value.expand())

    # g_id has only the differential specified
    g_1d = s.fxn('g_1d', None, g_1s.value.subs(f_1.dict).diff(x) * d(x) + g_1s.value.subs(f_1.dict).diff(y) * d(y))
    g_2d = s.fxn('g_2d', None, g_2s.value.subs(f_2.dict).diff(x) * d(x) + g_2s.value.subs(f_2.dict).diff(y) * d(y))
    g_3d = s.fxn('g_3d', None, g_3s.value.subs(f_3.dict).diff(x) * d(x) + g_3s.value.subs(f_3.dict).diff(y) * d(y))
    g_4d = s.fxn('g_4d', None, g_4s.value.subs(f_4.dict).diff(x) * d(x) + g_4s.value.subs(f_4.dict).diff(y) * d(y))

    test_passed = True
    
    test_passed = test_passed and bool(f_1.diff(x) == f_1.value.diff(x))
    test_passed = test_passed and bool(f_2.diff(x) == f_2.value.diff(x))
    test_passed = test_passed and bool(f_3.diff(x) == f_3.value.diff(x))
    test_passed = test_passed and bool(f_4.diff(x) == f_4.value.diff(x))
    test_passed = test_passed and bool(f_1.diff(y) == f_1.value.diff(y))
    test_passed = test_passed and bool(f_2.diff(y) == f_2.value.diff(y))
    test_passed = test_passed and bool(f_3.diff(y) == f_3.value.diff(y))
    test_passed = test_passed and bool(f_4.diff(y) == f_4.value.diff(y))

    test_passed = test_passed and bool(g_1e.diff(x).subs(f_1.dict) == g_1e.value.subs(f_1.dict).diff(x))
    test_passed = test_passed and bool(g_2e.diff(x).subs(f_2.dict) == g_2e.value.subs(f_2.dict).diff(x))
    test_passed = test_passed and bool(g_3e.diff(x).subs(f_3.dict) == g_3e.value.subs(f_3.dict).diff(x))
    test_passed = test_passed and bool(g_4e.diff(x).subs(f_4.dict) == g_4e.value.subs(f_4.dict).diff(x))
    test_passed = test_passed and bool(g_1e.diff(y).subs(f_1.dict) == g_1e.value.subs(f_1.dict).diff(y))
    test_passed = test_passed and bool(g_2e.diff(y).subs(f_2.dict) == g_2e.value.subs(f_2.dict).diff(y))
    test_passed = test_passed and bool(g_3e.diff(y).subs(f_3.dict) == g_3e.value.subs(f_3.dict).diff(y))
    test_passed = test_passed and bool(g_4e.diff(y).subs(f_4.dict) == g_4e.value.subs(f_4.dict).diff(y))

    test_passed = test_passed and bool(g_1s.diff(x).subs(f_1.dict) == g_1s.value.subs(f_1.dict).diff(x))
    test_passed = test_passed and bool(g_2s.diff(x).subs(f_2.dict) == g_2s.value.subs(f_2.dict).diff(x))
    test_passed = test_passed and bool(g_3s.diff(x).subs(f_3.dict) == g_3s.value.subs(f_3.dict).diff(x))
    test_passed = test_passed and bool(g_4s.diff(x).subs(f_4.dict) == g_4s.value.subs(f_4.dict).diff(x))
    test_passed = test_passed and bool(g_1s.diff(y).subs(f_1.dict) == g_1s.value.subs(f_1.dict).diff(y))
    test_passed = test_passed and bool(g_2s.diff(y).subs(f_2.dict) == g_2s.value.subs(f_2.dict).diff(y))
    test_passed = test_passed and bool(g_3s.diff(y).subs(f_3.dict) == g_3s.value.subs(f_3.dict).diff(y))
    test_passed = test_passed and bool(g_4s.diff(y).subs(f_4.dict) == g_4s.value.subs(f_4.dict).diff(y))

    test_passed = test_passed and bool(g_1e.diff(x).subs(f_1.dict) == g_1s.diff(x).subs(f_1.dict))
    test_passed = test_passed and bool(g_2e.diff(x).subs(f_2.dict) == g_2s.diff(x).subs(f_2.dict))
    test_passed = test_passed and bool(g_3e.diff(x).subs(f_3.dict) == g_3s.diff(x).subs(f_3.dict))
    test_passed = test_passed and bool(g_4e.diff(x).subs(f_4.dict) == g_4s.diff(x).subs(f_4.dict))
    test_passed = test_passed and bool(g_1e.diff(y).subs(f_1.dict) == g_1s.diff(y).subs(f_1.dict))
    test_passed = test_passed and bool(g_2e.diff(y).subs(f_2.dict) == g_2s.diff(y).subs(f_2.dict))
    test_passed = test_passed and bool(g_3e.diff(y).subs(f_3.dict) == g_3s.diff(y).subs(f_3.dict))
    test_passed = test_passed and bool(g_4e.diff(y).subs(f_4.dict) == g_4s.diff(y).subs(f_4.dict))

    test_passed = test_passed and bool(g_1d.diff(x).subs(f_1.dict) == g_1s.diff(x).subs(f_1.dict))
    test_passed = test_passed and bool(g_2d.diff(x).subs(f_2.dict) == g_2s.diff(x).subs(f_2.dict))
    test_passed = test_passed and bool(g_3d.diff(x).subs(f_3.dict) == g_3s.diff(x).subs(f_3.dict))
    test_passed = test_passed and bool(g_4d.diff(x).subs(f_4.dict) == g_4s.diff(x).subs(f_4.dict))
    test_passed = test_passed and bool(g_1d.diff(y).subs(f_1.dict) == g_1s.diff(y).subs(f_1.dict))
    test_passed = test_passed and bool(g_2d.diff(y).subs(f_2.dict) == g_2s.diff(y).subs(f_2.dict))
    test_passed = test_passed and bool(g_3d.diff(y).subs(f_3.dict) == g_3s.diff(y).subs(f_3.dict))
    test_passed = test_passed and bool(g_4d.diff(y).subs(f_4.dict) == g_4s.diff(y).subs(f_4.dict))
    
    return test_passed

def test_is_constant():
    s = PoblanoRing(['x', 'y'], {'x': '\\mathrm{x}'})
    x, y = s.vars
    d = s.d

    a = SageVariable('a')
    b = SageVariable('b')

    af = s.fxn('a_f')
    bf = s.fxn('b_f')

    c = s.fxn('c', 2 * sqrt(2) * af + bf)

    # we test across constants that are numbers (2, 3), numeric symbols (sqrt(2), sqrt(3)), Sage variables (a, b), and constant PoblanoExpression objects (af, bf)
    f_1 = s.fxn('f_1', 2 * x + 3 * y)
    f_2 = s.fxn('f_2', sqrt(2) * x + sqrt(3) * y)
    f_3 = s.fxn('f_3', a * x + b * y)
    f_4 = s.fxn('f_4', af * x + bf * y)

    # g_is is given in a condensed form
    g_1s = s.fxn('g_1s', (x * y * (2 * sqrt(2) * c * f_1 + f_1 ** 2 / y) + exp(f_1 * x)))
    g_2s = s.fxn('g_2s', (x * y * (2 * sqrt(2) * c * f_2 + f_2 ** 2 / y) + exp(f_2 * x)))
    g_3s = s.fxn('g_3s', (x * y * (2 * sqrt(2) * c * f_3 + f_3 ** 2 / y) + exp(f_3 * x)))
    g_4s = s.fxn('g_4s', (x * y * (2 * sqrt(2) * c * f_4 + f_4 ** 2 / y) + exp(f_4 * x)))

    # g_ie is expanded
    g_1e = s.fxn('g_1e', g_1s.value.expand())
    g_2e = s.fxn('g_2e', g_2s.value.expand())
    g_3e = s.fxn('g_3e', g_3s.value.expand())
    g_4e = s.fxn('g_4e', g_4s.value.expand())

    # g_id has only the differential specified
    g_1d = s.fxn('g_1d', None, g_1s.value.subs(f_1.dict).diff(x) * d(x) + g_1s.value.subs(f_1.dict).diff(y) * d(y))
    g_2d = s.fxn('g_2d', None, g_2s.value.subs(f_2.dict).diff(x) * d(x) + g_2s.value.subs(f_2.dict).diff(y) * d(y))
    g_3d = s.fxn('g_3d', None, g_3s.value.subs(f_3.dict).diff(x) * d(x) + g_3s.value.subs(f_3.dict).diff(y) * d(y))
    g_4d = s.fxn('g_4d', None, g_4s.value.subs(f_4.dict).diff(x) * d(x) + g_4s.value.subs(f_4.dict).diff(y) * d(y))
    
    test_passed = True
    
    test_passed = test_passed and bool(s.is_constant(2) is True)
    test_passed = test_passed and bool(s.is_constant(3) is True)
    test_passed = test_passed and bool(s.is_constant(sqrt(2)) is True)
    test_passed = test_passed and bool(s.is_constant(sqrt(3)) is True)
    test_passed = test_passed and bool(s.is_constant(af) is True)
    test_passed = test_passed and bool(s.is_constant(bf) is True)
    test_passed = test_passed and bool(s.is_constant(c) is True)
    test_passed = test_passed and bool(s.is_constant(x) is False)
    test_passed = test_passed and bool(s.is_constant(y) is False)
    test_passed = test_passed and bool(s.is_constant(f_1) is False)
    test_passed = test_passed and bool(s.is_constant(f_2) is False)
    test_passed = test_passed and bool(s.is_constant(f_3) is False)
    test_passed = test_passed and bool(s.is_constant(f_4) is False)
    test_passed = test_passed and bool(s.is_constant(g_1e) is False)
    test_passed = test_passed and bool(s.is_constant(g_2e) is False)
    test_passed = test_passed and bool(s.is_constant(g_3e) is False)
    test_passed = test_passed and bool(s.is_constant(g_4e) is False)
    test_passed = test_passed and bool(s.is_constant(g_1s) is False)
    test_passed = test_passed and bool(s.is_constant(g_2s) is False)
    test_passed = test_passed and bool(s.is_constant(g_3s) is False)
    test_passed = test_passed and bool(s.is_constant(g_4s) is False)
    test_passed = test_passed and bool(s.is_constant(g_1d) is False)
    test_passed = test_passed and bool(s.is_constant(g_2d) is False)
    test_passed = test_passed and bool(s.is_constant(g_3d) is False)
    test_passed = test_passed and bool(s.is_constant(g_4d) is False)
    
    return test_passed

if test_constant_poblano_expression_diff_and_jump() and \
       test_poblano_expression_construction_options() and \
       test_poblano_expression_vs_space_diff_and_jump() and \
       test_jump_expansion_accuracy_simple_poblano_expression() and \
       test_jump_expansion_accuracy_poblano_expression_and_constants() and \
       test_jump_expansion_accuracy_binary_operations() and \
       test_jump_expansion_accuracy_integer_powers() and \
       test_jump_expansion_accuracy_noninteger_powers() and \
       test_jump_expansion_accuracy_unary_operations() and \
       test_jump_expansion_accuracy_polynomial_consistency() and \
       test_differentiation_over_constants_and_expands() and \
       test_is_constant():
    print('all tests passed!')
    exit(0)
else:
    print('not all tests passed!')
    exit(1)