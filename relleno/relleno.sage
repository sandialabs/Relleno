# Relleno - a SageMath library for automatic jump expansion and convenient calculus                        
# Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).                       
#                                                                                                        
# This program is free software: you can redistribute it and/or modify                                     
# it under the terms of the GNU General Public License as published by                                     
# the Free Software Foundation, either version 3 of the License, or                                        
# (at your option) any later version.                                                                      
#                                                                                                          
# This program is distributed in the hope that it will be useful,                                          
# but without any warranty; without even the implied warranty of                                           
# merchantability or fitness for a particular purpose.  See the                                            
# GNU General Public License for more details.                                                             
#                                                                                                          
# You should have received a copy of the GNU General Public License                                        
# along with this program.  If not, see <http://www.gnu.org/licenses/>.                                    
#                                                                                                          
# Questions? Contact Mike Hansen (mahanse@sandia.gov)      

from sage.all import show, assume, diff, sqrt, exp, latex, log, prod, vector, matrix, identity_matrix
from sage.all import var as SageVariable
from sage.all import SR as Symbols, RR as Reals, ZZ as Integers
from sage.all import SageObject


class RellenoExpression(SageObject):
    """The class that allows `hidden' variable values, derivatives, and jump expansions.

    RellenoExpression objects can be used in most Sage expressions, although occasionally (e.g., exp, log, cos, etc.)
    you'll need to grab the symbol off of the object (e.g., exp(a.symbol)).

    Note that it is better to build RellenoExpression through the fxn() method of a RellenoRing

    **Constructor** : specify a relleno ring, the symbol, and optional value, total differential, jump dictionary, and LaTeX representation

    :param ring: the RellenoRing on which this expression exists
    :param symbol: the symbol (a string) to use in SageMath expressions
    :param value: the value of this expression (a sage.symbolic.expression.Expression, optional, default: None)
    :param total_diff: the total differential of this expression (a sage.symbolic.expression.Expression, optional, default: None)
    :param jump_dict: the dictionary mapping variables to jump ratios (optional, default: None)
    :param latex_name: the LaTeX representation of this expression for Sage's pretty print (show) typesetting

    Useful combinations of constructor parameters (all combinations include the ring and optionally the latex_name):

    1. symbol
        With only the symbol provided, the RellenoExpression is treated as a constant - effectively just makes a Sage variable
    2. symbol, value
        In this case the total differential and jump expansion are computed on construction from the given value
    3. symbol, total_diff
        Here the jump expansion is treated as zero, as it cannot be inferred from the derivative alone
    4. symbol, jump_dict
        With the jump expansion given, the derivative is inferred from the jump expansion at the limit of zero jumps
    5. symbol, value, total_diff
        Here the jump dictionary is computed from the given value, there is no check that the specified differential matches the computed jump expansion at the zero jump limit
    6. symbol, value, jump_dict
        Here the total differential is computed from the given value
    7. symbol, value, total_diff, jump_dict
        Here everything is retained in its given form

    """
    def __init__(self, ring, symbol, value=None, total_diff=None, jump_dict=None, latex_name=None):
        self._ring = ring

        if latex_name is None:
            self._symbol = SageVariable(symbol)
        else:
            self._symbol = SageVariable(symbol, latex_name=latex_name)

        self._value = value if value is not None else self._symbol
        self._dict = dict({self._symbol: self._value})

        if total_diff is not None:
            self._total_diff = total_diff
        else:
            self._total_diff = 0
            for v in self._ring.vars:
                self._total_diff += self._ring.diff(self._value, v) * self._ring.d(v)

        if jump_dict is None:
            self._jump_dict = dict()
            if value is None:
                for v in self._ring.vars:
                    self._jump_dict[v] = 0
            else:
                jump_expr = self._ring._build_jump_expr(value)
                for v in self._ring.vars:
                    self._jump_dict[v] = jump_expr.jump_coeff(v)
        else:
            self._jump_dict = dict()
            for v in self._ring.vars:
                if v not in jump_dict:
                    self._jump_dict[v] = 0
                else:
                    self._jump_dict[v] = jump_dict[v]

        if jump_dict is not None and total_diff is None and value is None:
            for v in self._ring.vars:
                self._total_diff += self._ring.subs_for_consistency(self._jump_dict[v]) * self._ring.d(v)

        self._has_zero_diff = bool(self._total_diff == 0)
        self._has_zero_jump = bool(sum([self._jump_dict[v] for v in self._ring.vars]) == 0)
        self._is_constant = self._has_zero_diff and self._has_zero_jump

    def diff(self, v):
        """Compute the derivative of this expression with respect to a variable of the RellenoRing"""
        return self._total_diff.diff(self._ring.d(v))

    @property
    def symbol(self):
        """Obtain the symbol of this expression"""
        return self._symbol

    @property
    def value(self):
        """Obtain this expression's value"""
        return self._value

    @property
    def latex_name(self):
        """Obtain the LaTeX representation of this expression"""
        return self._latex_name

    def is_constant(self):
        """Check whether or not this expression is a constant"""
        return self._is_constant

    @property
    def total_diff(self):
        """Obtain the total differential of this expression"""
        return self._total_diff

    @property
    def jump_dict(self):
        """Obtain the jump expansion (dictionary) of this expression"""
        return self._jump_dict

    @property
    def dict(self):
        """Obtain the symbol:value dictionary of this expression (useful for substitutions)"""
        return self._dict

    def jump_coeff(self, v):
        """Obtain the jump coefficient (ratio) of this expression with respect to a given variable of the RellenoRing"""
        return self._jump_dict[v]

    def has_zero_jump(self, ):
        """Check whether or not this expression has a trivial jump expansion"""
        return self._has_zero_jump

    def __eq__(self, other):
        if isinstance(other, RellenoExpression):
            return self._dict == other.dict and self._total_diff == other.total_diff
        else:
            return False

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'RellenoExpr: (' + latex(self.symbol) + ': ' + \
               latex(self.value) + ' (' + ('constant))' if self.is_constant() else 'nonconstant))')

    def __hash__(self):
        return str(self)

    def _latex_(self):
        return latex(self._symbol)

    def __neg__(self):
        return -self.symbol

    def __add__(self, other):
        return self._symbol + other

    __radd__ = __add__

    def __mul__(self, other):
        return self._symbol * other

    __rmul__ = __mul__

    def __pow__(self, other):
        return self._symbol ** other

    def __rpow__(self, other):
        return other ** self._symbol

    def __sub__(self, other):
        return self._symbol - other

    def __rsub__(self, other):
        return other - self._symbol

    def __div__(self, other):
        return self._symbol / other

    def __truediv__(self, other):
        return self._symbol / other

    def __rdiv__(self, other):
        return other / self._symbol

    def __rtruediv__(self, other):
        return other / self._symbol


class RellenoRing(SageObject):
    """The class that facilitates differentiation, averaging, and jump expansion of RellenoExpression and Sage expressions

    **Constructor**: build a RellenoRing from a set of fundamental variables

    :param name: the name of this RellenoRing object, must be a valid Python identifier (no spaces or dashes)
    :param base_variables: a list of strings representing the fundamental variables of this ring
    :param latex_name_dict: a dictionary that allows special LaTeX representations of the variables (optional, default: None)
    """
    def __init__(self, name, base_variables, latex_name_dict=None):
        self._name = name
        if latex_name_dict is None:
            self._vars = [SageVariable(b) for b in base_variables]
        else:
            self._vars = [SageVariable(b, latex_name=latex_name_dict[b] if b in latex_name_dict else b) for b in
                          base_variables]
        self._vdv = {v: SageVariable('d' + str(v), latex_name='\\mathrm{d}' + str(v)) for v in self._vars}
        self._fxns = dict()
        self._tmps = dict()
        self._invjumptmps = dict()

        self._mark_dict = {'arithmetic': '\\overline', 'logarithmic': '\\widehat', 'inverse': '\\widetilde'}
        self._tmpsymbol = 'xi_'

        self._oper_dict = dict()
        operator_template_1 = SageVariable('operator_template_1')
        operator_template_2 = SageVariable('operator_template_2')
        self._oper_dict[log(operator_template_1).operator()] = 'log'
        self._oper_dict[(operator_template_1 * operator_template_2).operator()] = 'mlt'
        self._oper_dict[(operator_template_1 + operator_template_2).operator()] = 'add'
        self._oper_dict[(operator_template_1 ** operator_template_2).operator()] = 'pow'

    @property
    def vars(self):
        """Obtain the fundamental variables that define this ring"""
        return self._vars

    def d(self, v):
        """Obtain the differential of a given variable, e.g., 'dx' or 'dy'"""
        return self._vdv[v]

    @property
    def temporaries(self):
        """Obtain the temporary variables defined in the process of computing jump expansions, averages, etc."""
        return self._tmps

    def fxn(self, symbol, value=None, total_diff=None, jump_dict=None, latex_name=None):
        """Build a new function (or constant) on this ring

        :param symbol:the symbol (a string) to use in SageMath expressions
        :param value: the value of this expression (a sage.symbolic.expression.Expression, optional, default: None)
        :param total_diff: the total differential of this expression (a sage.symbolic.expression.Expression, optional, default: None)
        :param jump_dict: the dictionary mapping variables to jump ratios (optional, default: None)
        :param latex_name: the LaTeX representation of this expression for Sage's pretty print (show) typesetting
        :return: a RellenoExpression

        Useful combinations of constructor parameters (all combinations include the ring and optionally the latex_name):

        1. symbol
            With only the symbol provided, the RellenoExpression is treated as a constant - effectively just makes a Sage variable
        2. symbol, value
            In this case the total differential and jump expansion are computed on construction from the given value
        3. symbol, total_diff
            Here the jump expansion is treated as zero, as it cannot be inferred from the derivative alone
        4. symbol, jump_dict
            With the jump expansion given, the derivative is inferred from the jump expansion at the limit of zero jumps
        5. symbol, value, total_diff
            Here the jump dictionary is computed from the given value, there is no check that the specified differential matches the computed jump expansion at the zero jump limit
        6. symbol, value, jump_dict
            Here the total differential is computed from the given value
        7. symbol, value, total_diff, jump_dict
            Here everything is retained in its given form
        """
        new_expr = RellenoExpression(self, symbol, value, total_diff, jump_dict, latex_name)
        self._fxns[new_expr.symbol] = new_expr
        return self._fxns[new_expr.symbol]

    def _get_all_symbolic_operands(self, an_expr):
        if isinstance(an_expr, RellenoExpression):
            return set({an_expr.symbol})
        elif an_expr in self._fxns:
            return set({an_expr})
        else:
            if not len(an_expr.operands()):
                return set({an_expr})
            else:
                full_list = set()
                ops = an_expr.operands()
                for op in ops:
                    if not len(op.operands()):
                        if not op.is_numeric():
                            full_list.add(op)
                    else:
                        full_list.update(self._get_all_symbolic_operands(op))
                return full_list

    def diff(self, expr, prim_var):
        """Compute the derivative of a Sage expression with respect to a fundamental variable of this ring"""
        dexpr = 0
        if expr.parent() == Symbols or isinstance(expr, RellenoExpression):
            all_ops = self._get_all_symbolic_operands(expr)
            for op in all_ops:
                if op in self._fxns:
                    if isinstance(expr, RellenoExpression):
                        dexpr += self._fxns[op].diff(prim_var)
                    else:
                        dexpr += expr.diff(op) * self._fxns[op].diff(prim_var)
                else:
                    if expr.is_numeric():
                        dexpr += 0
                    else:
                        dexpr += expr.diff(op) * op.diff(prim_var)
        return dexpr

    def is_constant(self, value):
        """Check whether or not a Sage expression is constant"""
        if value in Reals or value in Integers or (value in Symbols and value.is_numeric()):
            return True
        else:
            if isinstance(value, RellenoExpression):
                return value.is_constant()
            elif value in self._fxns:
                return self._fxns[value].is_constant()
            else:

                return bool(sum([self.diff(value, v) * self.d(v) for v in self._vars]) == 0) and self.has_zero_jump(
                    value)

    def avg(self, expr, type='arithmetic'):
        """ Obtain the average of a Sage expression

        :param expr: the Sage expression
        :param type: the type of average (usually 'arithmetic', 'inverse', 'logarithmic', or an operator)
        :return: either a RellenoExpression or a Sage expression (depending upon the input)

        This function will aggressively try to remove constants from arithmetic averages as much as possible.
        For instance, avg(2f) is returned as 2avg(f)
        """
        if self.is_constant(expr):
            return expr
        else:
            for t in self._tmps:
                if bool(expr == self._tmps[t]['expression']) and type == self._tmps[t]['avg type']:
                    return self._tmps[t]['temporary']

            var_name = self._name + '_' + self._tmpsymbol + str(len(self._tmps))

            latex_expr = latex(expr)
            if isinstance(expr, RellenoExpression):
                latex_expr = latex(expr.symbol)
            if type in self._mark_dict:
                latex_name = self._mark_dict[type] + '{' + latex_expr + '}\,'
            else:
                latex_name = '{\\left[\\left[' + latex_expr + '\\right]\\right]_{\\mathrm{' + type + '}}}\,'

            if isinstance(expr, RellenoExpression):
                if var_name not in self._tmps:
                    var_val = SageVariable(var_name + '_value', latex_name=latex_name)
                    self._tmps[var_name] = {'temporary': self.fxn(var_name, var_val, latex_name=latex_name),
                                            'expression': expr,
                                            'avg type': type}
                return self._tmps[var_name]['temporary']
            elif expr in self._fxns:
                return self.avg(self._fxns[expr], type)
            else:
                operands, operator = expr.operands(), expr.operator()
                if type == 'arithmetic' and (
                        operator in self._oper_dict and self._oper_dict[operator] in ['mlt', 'add']):
                    cons = []
                    nons = []
                    for op in operands:
                        if self.is_constant(op):
                            cons.append(op)
                        else:
                            nons.append(op)
                    if not cons:
                        if var_name not in self._tmps:
                            var_val = SageVariable(var_name + '_value', latex_name=latex_name)
                            self._tmps[var_name] = {'temporary': self.fxn(var_name, var_val, latex_name=latex_name),
                                                    'expression': expr,
                                                    'avg type': type}
                        return self._tmps[var_name]['temporary']
                    else:
                        if not nons:
                            return expr
                        else:
                            if self._oper_dict[operator] == 'mlt':
                                return prod(cons) * self.avg(prod(nons), type)
                            elif self._oper_dict[operator] == 'add':
                                return sum(cons) + sum([self.avg(non, type) for non in nons])

                else:
                    if var_name not in self._tmps:
                        var_val = SageVariable(var_name + '_value', latex_name=latex_name)
                        self._tmps[var_name] = {'temporary': self.fxn(var_name, var_val, latex_name=latex_name),
                                                'expression': expr,
                                                'avg type': type}
                    return self._tmps[var_name]['temporary']

    def subs_for_consistency(self, expr):
        """Substitute temporaries from averages and jump calculations and reduce a two-point expression into its limit of zero jump."""
        subs_dict = dict()
        for t in self.temporaries:
            temp_expr = self.temporaries[t]['expression']
            temp_type = self.temporaries[t]['avg type']
            temp = self.temporaries[t]['temporary'].symbol
            if temp_type == 'arithmetic':
                if isinstance(temp_expr, RellenoExpression):
                    subs_dict[temp] = temp_expr.symbol
                else:
                    subs_dict[temp] = temp_expr
            elif temp_type == 'inverse':
                subs_dict[temp] = -temp_expr * temp_expr
            elif temp_type == 'logarithmic':
                if isinstance(temp_expr, RellenoExpression):
                    subs_dict[temp] = temp_expr.symbol
                else:
                    subs_dict[temp] = temp_expr
            else:
                bleh = SageVariable('bleh')
                operator_derivative = eval(temp_type + '(' + str(bleh) + ').diff(' + str(bleh) + ')')
                if isinstance(expr, RellenoExpression):
                    expr_to_subs = temp_expr.symbol
                else:
                    expr_to_subs = temp_expr
                subs_dict[temp] = operator_derivative.subs({bleh: expr_to_subs})
        if isinstance(expr, RellenoExpression):
            subbed_expr = expr.value
        elif expr in self._fxns:
            subbed_expr = self._fxns[expr].value
        else:
            subbed_expr = expr
        subbed_expr = subbed_expr.subs(subs_dict)
        for t in self._invjumptmps:
            subbed_expr = subbed_expr.subs(self._invjumptmps[t].dict)
        return subbed_expr

    def jump_coeff(self, expr, jump_var):
        """Compute the jump coefficient (ratio) of a Sage expression with respect to a fundamental variable of this ring"""
        return self._build_jump_expr(expr).jump_coeff(jump_var)

    def has_zero_jump(self, expr):
        """Check whether or not a Sage expression has a trivial jump expansion"""
        return self._build_jump_expr(expr).has_zero_jump()

    def _build_jump_expr(self, expr):
        if expr in Reals or expr in Integers or (expr in Symbols and expr.is_numeric()):
            return _ConstantExpr(expr, self)
        if isinstance(expr, RellenoExpression):
            return _SpecialJumpExpr(expr, self)
        elif expr in self._fxns:
            return _SpecialJumpExpr(self._fxns[expr], self)
        else:
            operands = expr.operands()
            operator = expr.operator()

            def throw_unknown(opstr):
                raise ValueError('unknown operation type detected: ' + opstr + ' and general form not implemented yet')

            if operator is None:
                if expr in self._vars:
                    return _VariableExpr(expr, self)
                elif expr in self._fxns:
                    return _SpecialJumpExpr(self._fxns[expr], self)
                else:
                    return _ConstantExpr(expr, self)
            else:
                if operator not in self._oper_dict:
                    return _GeneralExpr(operator, self._build_jump_expr(operands[0]), self)
                else:
                    opstr = self._oper_dict[operator]
                    if opstr == 'add':
                        expressions = []
                        for operand in operands:
                            expressions.append(self._build_jump_expr(operand))
                        return _SumExpr(expressions, self)
                    elif opstr == 'mlt':
                        expressions = []
                        for operand in operands:
                            expressions.append(self._build_jump_expr(operand))
                        return _ProdExpr(expressions, self)
                    elif opstr == 'log':
                        return _LogExpr(self._build_jump_expr(operands[0]), self)
                    elif opstr == 'pow':
                        quantity = operands[0]
                        exponent = operands[1]
                        if exponent > 0:
                            if exponent in Integers:
                                if not int(exponent) % 2:
                                    return _SquareExpr(self._build_jump_expr(quantity ** (exponent / 2)), self)
                                else:
                                    return _ProdExpr([self._build_jump_expr(quantity),
                                                      _SquareExpr(self._build_jump_expr(quantity ** ((exponent - 1) / 2)),
                                                                  self)], self)
                                return _ProdExpr(expressions, self)
                            else:
                                if exponent == 1 / 2:
                                    return _SqrtExpr(self._build_jump_expr(quantity), self)
                                else:
                                    return self._build_jump_expr(exp(exponent * log(quantity)))
                        elif exponent < 0:
                            if exponent in Integers:
                                if exponent == -1:
                                    return _InvExpr(self._build_jump_expr(quantity), self)
                                else:
                                    invq = self.fxn('inverse_' + str(quantity), 1 / quantity,
                                                    latex_name='\\left(1/' + latex(quantity) + '\\right)')
                                    if isinstance(quantity, RellenoExpression):
                                        self._invjumptmps[quantity.symbol] = invq
                                    else:
                                        self._invjumptmps[quantity] = invq
                                    n = -exponent
                                    return self._build_jump_expr(invq ** n)
                            else:
                                if exponent == -1 / 2:
                                    return _InvExpr(self._build_jump_expr(quantity ** (-exponent)), self)
                                else:
                                    return self._build_jump_expr(exp(exponent * log(quantity)))
                        else:
                            return self._build_jump_expr(exp(exponent * log(quantity)))


class _JumpExpr(object):
    def __init__(self, sage_expr, relleno_ring):
        self._sage_expr = sage_expr
        self._relleno_ring = relleno_ring

    def is_constant(self):
        return False

    def is_variable(self):
        return False

    @property
    def sage_expr(self):
        return self._sage_expr

    @property
    def relleno_ring(self):
        return self._relleno_ring

    def has_zero_jump(self):
        return False

    def jump_coeff(self, variable):
        raise Exception('Called method jump_coeff(variable) from base class (JumpExpr)!')


class _ConstantExpr(_JumpExpr):
    def __init__(self, sage_expr, relleno_ring):
        _JumpExpr.__init__(self, sage_expr, relleno_ring)

    def is_constant(self):
        return True

    def has_zero_jump(self):
        return True

    def jump_coeff(self, variable):
        return 0


class _VariableExpr(_JumpExpr):
    def __init__(self, sage_expr, relleno_ring):
        _JumpExpr.__init__(self, sage_expr, relleno_ring)

    def is_variable(self):
        return True

    def has_zero_jump(self):
        return False

    def jump_coeff(self, variable):
        return 1 if variable == self.sage_expr else 0


class _SpecialJumpExpr(_JumpExpr):
    def __init__(self, relleno_expr, relleno_ring):
        self._relleno_expr = relleno_expr
        _JumpExpr.__init__(self, relleno_expr.symbol, relleno_ring)

    def jump_coeff(self, variable):
        return self._relleno_expr.jump_coeff(variable)

    def has_zero_jump(self):
        return self._relleno_expr.has_zero_jump()


class _SumExpr(_JumpExpr):
    def __init__(self, tree_operands, relleno_ring):
        summation = tree_operands[0].sage_expr
        for to in tree_operands[1:]:
            summation += to.sage_expr
        _JumpExpr.__init__(self, summation, relleno_ring)
        self.tree_operands = tree_operands

    def jump_coeff(self, variable):
        return sum([top.jump_coeff(variable) for top in self.tree_operands])

    def has_zero_jump(self):
        return bool(all([top.has_zero_jump() for top in self.tree_operands]))


class _ProdExpr(_JumpExpr):
    def __init__(self, tree_operands, relleno_ring):
        product = tree_operands[0].sage_expr
        for to in tree_operands[1:]:
            product *= to.sage_expr
        _JumpExpr.__init__(self, product, relleno_ring)
        self.tree_operands = tree_operands
        self.is_binary = len(self.tree_operands) == 2

    def jump_coeff(self, variable):
        if self.is_binary:
            avg0 = self._relleno_ring.avg(self.tree_operands[0].sage_expr, 'arithmetic')
            avg1 = self._relleno_ring.avg(self.tree_operands[1].sage_expr, 'arithmetic')
            return avg0 * self.tree_operands[1].jump_coeff(variable) + avg1 * self.tree_operands[0].jump_coeff(variable)
        else:
            # split into binary (left->right linked list) and then reapply
            left_expr = self.tree_operands[0]
            right_expr = _ProdExpr(self.tree_operands[1:], self.relleno_ring)
            avg0 = self._relleno_ring.avg(left_expr.sage_expr, 'arithmetic')
            avg1 = self._relleno_ring.avg(right_expr.sage_expr, 'arithmetic')
            return avg0 * right_expr.jump_coeff(variable) + avg1 * left_expr.jump_coeff(variable)

    def has_zero_jump(self):
        return bool(all([top.has_zero_jump() for top in self.tree_operands]))


class _LogExpr(_JumpExpr):
    def __init__(self, tree_operand, relleno_ring):
        _JumpExpr.__init__(self, log(tree_operand.sage_expr), relleno_ring)
        self.tree_operand = tree_operand

    def jump_coeff(self, variable):
        if self.tree_operand.is_constant():
            return 0
        elif self.tree_operand.is_variable():
            if self.tree_operand.sage_expr == variable:
                log_avg = self._relleno_ring.avg(self.tree_operand.sage_expr, 'logarithmic')
                return 1 / log_avg
            else:
                return 0
        else:
            log_avg = self._relleno_ring.avg(self.tree_operand.sage_expr, 'logarithmic')
            return 1 / log_avg * self.tree_operand.jump_coeff(variable)

    def has_zero_jump(self):
        return bool(self.tree_operand.has_zero_jump())


class _InvExpr(_JumpExpr):
    def __init__(self, tree_operand, relleno_ring):
        _JumpExpr.__init__(self, 1 / tree_operand.sage_expr, relleno_ring)
        self.tree_operand = tree_operand

    def jump_coeff(self, variable):
        if self.tree_operand.is_constant():
            return 0
        elif self.tree_operand.is_variable():
            if self.tree_operand.sage_expr == variable:
                invavg = self._relleno_ring.avg(self.tree_operand.sage_expr, 'inverse')
                return 1 / invavg
            else:
                return 0
        else:
            invavg = self._relleno_ring.avg(self.tree_operand.sage_expr, 'inverse')
            return self.tree_operand.jump_coeff(variable) / invavg

    def has_zero_jump(self):
        return bool(self.tree_operand.has_zero_jump())


class _GeneralExpr(_JumpExpr):
    def __init__(self, operator, tree_operand, relleno_ring):
        _JumpExpr.__init__(self, operator(tree_operand.sage_expr), relleno_ring)
        self.tree_operand = tree_operand
        self.operator = operator

    def jump_coeff(self, variable):
        if self.tree_operand.is_constant():
            return 0
        elif self.tree_operand.is_variable():
            if self.tree_operand.sage_expr == variable:
                genavg = self._relleno_ring.avg(self.tree_operand.sage_expr, str(self.operator))
                return 1 / genavg
            else:
                return 0
        else:
            genavg = self._relleno_ring.avg(self.tree_operand.sage_expr, str(self.operator))
            return self.tree_operand.jump_coeff(variable) / genavg

    def has_zero_jump(self):
        return bool(self.tree_operand.has_zero_jump())


class _SqrtExpr(_JumpExpr):
    def __init__(self, tree_operand, relleno_ring):
        _JumpExpr.__init__(self, sqrt(tree_operand.sage_expr), relleno_ring)
        self.tree_operand = tree_operand

    def jump_coeff(self, variable):
        if self.tree_operand.is_constant():
            return 0
        elif self.tree_operand.is_variable():
            if self.tree_operand.sage_expr == variable:
                avgsqrt = self._relleno_ring.avg(sqrt(self.tree_operand.sage_expr), 'arithmetic')
                return 1 / (2 * avgsqrt)
            else:
                return 0
        else:
            avgsqrt = self._relleno_ring.avg(sqrt(self.tree_operand.sage_expr), 'arithmetic')
            return self.tree_operand.jump_coeff(variable) / (2 * avgsqrt)

    def has_zero_jump(self):
        return bool(self.tree_operand.has_zero_jump())


class _SquareExpr(_JumpExpr):
    def __init__(self, tree_operand, relleno_ring):
        _JumpExpr.__init__(self, tree_operand.sage_expr ** 2, relleno_ring)
        self.tree_operand = tree_operand

    def jump_coeff(self, variable):
        if self.tree_operand.is_constant():
            return 0
        elif self.tree_operand.is_variable():
            if self.tree_operand.sage_expr == variable:
                avg = self._relleno_ring.avg(self.tree_operand.sage_expr, 'arithmetic')
                return 2 * avg
            else:
                return 0
        else:
            avg = self._relleno_ring.avg(self.tree_operand.sage_expr, 'arithmetic')
            return 2 * avg * self.tree_operand.jump_coeff(variable)

    def has_zero_jump(self):
        return bool(self.tree_operand.has_zero_jump())