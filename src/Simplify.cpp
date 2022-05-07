#include "Simplify.hpp"
#include "Traverse.hpp"
#include "Utilities.hpp"
#include "Hashing.hpp"
#include "CodeGenerator.hpp"
#include <unordered_map>
#include <unordered_set>
#include <set>

using namespace std;

namespace Sym {

Symbolic makeMul(const vector<Symbolic>& childs) {
	if (childs.size() == 0) return Symbolic(1.); // not defined

	vector<Symbolic> childs2;
	double f = 1.;

	for (auto& c : childs) {
		if (c.op() == CONST) {
			f *= c.constant();
		}
		else if (isSmallConstant(c.ahash())) {
			f *= (double)c.ahash();
		}
		else childs2.push_back(c);
	}

	Symbolic ret;

	if (f == .0) return Symbolic(.0);
	if (f != 1.) childs2.push_back(Symbolic(f));
	if (childs2.size() == 0) return Symbolic(1.);
	if (childs2.size() == 1) ret = childs2.front();
	else ret = Symbolic(MUL, childs2);

	return ret;
}

Symbolic makeAddition(vector<Symbolic> operands) {

	double constant = .0;

	operands.erase(remove_if(operands.begin(), operands.end(), [&](const auto& x) {

		if (x.op() == CONST) {
			constant += x.constant();
			return true;
		}
		else return false;

		}), operands.end());

	if (constant != 0) operands.push_back(constant);

	if (operands.size() == 0) return .0;
	if (operands.size() == 1) return operands.front();

	return Symbolic(ADD, operands);
}

Symbolic makeSymbolic(const unsigned char op, vector<Symbolic> operands) {

	Symbolic ret;

	switch (op) {
	case ADD: ret = makeAddition(operands);
		break;

	case MUL:
	case MULC:
		ret = makeMul(operands);
		break;

	default: ret = Symbolic(op, operands);
	}

	if (isSmallConstant(ret.ahash())) return (double)ret.ahash();
	else return ret;
}

Symbolic flattenAdditions(const Symbolic& expr) {

	if (expr.op() != ADD) return expr;

	vector<Symbolic> operands;

	preOrderTraverse(expr, [&](const Symbolic& x) {

		if (!x.ahash()) return false;

		if (x.op() != ADD) {
			operands.push_back(x);
			return false;
		}

		return true;
		});

	auto ret = Symbolic(ADD, operands);
	if (ret.ahash() == .0) ret = .0;
	else if (ret.numChilds() == 1) ret = Symbolic(ret[0]);

	return ret;
}

Symbolic flattenMultiplications(const Symbolic& expr) {

	if (expr.op() != MUL && expr.op() != MULC) return expr;

	vector<Symbolic> operands;

	preOrderTraverse(expr, [&](const Symbolic& x) {

		if (isSmallConstant(x.ahash())) {
			if (x.ahash() != 1) operands.push_back((double)x.ahash());
			return false;
		}

		if (x.op() != MUL && x.op() != MULC) {
			operands.push_back(x);
			return false;
		}

		return true;
		});

	return makeSymbolic(MUL, operands);
}

Symbolic makePower(const Symbolic& x, const double d) {
	if (x.op() == CONST) return Symbolic(pow(x.constant(), d));
	else if (x.op() == POW && x[1].op() == CONST) return Symbolic(POW, x[0], Symbolic(d * x[1].constant()));
	else return Symbolic(POW, x, Symbolic(d));
}

bool isMul(const Symbolic& x) {
	return x.op() == MUL || x.op() == MULC;
}

Symbolic flattenDivAndMul(const Symbolic& expr) {

	vector<Symbolic> childs;

	switch (expr.op()) {
	case MUL:
	case MULC:

		for (auto& c : mapData(expr, [](const auto& x) {return flattenDivAndMul(x); })) {

			if (isMul(c)) for (auto& cc : c) childs.push_back(cc);
			else childs.push_back(c);
		}

		break;

	case SQRT: {
		auto rad = flattenDivAndMul(expr[0]);
		if (isMul(rad)) for (auto& c : rad) childs.push_back(makePower(c, 0.5));
		else childs.push_back(makePower(rad, 0.5));
		break;
	}

	case DIV: {
		auto nom = flattenDivAndMul(expr[0]);
		auto den = flattenDivAndMul(expr[1]);

		if (isMul(nom)) copy(nom.begin(), nom.end(), back_inserter(childs));
		else childs.push_back(nom);

		if (isMul(den)) for (auto& c : den) childs.push_back(makePower(c, -1.));
		else childs.push_back(makePower(den, -1.));
		break;
	}

	default: return expr;

	}

	return makeMul(childs);
}

Symbolic optimizeProducts(const Symbolic& x) {

	auto xf = flattenDivAndMul(x);
	if (!isMul(xf)) return x;

	map<Symbolic, double, AlgebraicHashFunctor> factors;

	for (auto& c : xf) {
		if (c.op() == POW) factors[c[0]] += c[1].constant();
		else factors[c] += 1.;
	}

	vector<Symbolic> nom, den;
	vector<Symbolic> nomR, denR;

	double f = 1.;

	for (auto& c : factors) {

		if (c.first.op() == CONST) {
			f *= pow(c.first.constant(), c.second);
		}
		else if (c.second < 0) {
			if (floor(c.second) == c.second) {
				for (int i = 0; i < floor(-c.second); ++i) den.push_back(c.first);
			}
			else if (floor(2 * c.second) == 2 * c.second) {
				for (int i = 0; i < floor(-c.second); ++i) den.push_back(c.first);
				denR.push_back(c.first);
			}
			else assert(0);

		}
		else if (c.second > 0) {
			if (floor(c.second) == c.second) {
				for (int i = 0; i < floor(c.second); ++i) nom.push_back(c.first);
			}
			else if (floor(2 * c.second) == 2 * c.second) {
				for (int i = 0; i < floor(c.second); ++i) nom.push_back(c.first);
				nomR.push_back(c.first);
			}
			else assert(0);
		}
	}

	Symbolic ret;

	if (!nomR.empty()) {
		if (denR.empty()) {
			nom.emplace_back(SQRT, makeMul(nomR));
		}
		else {
			nom.emplace_back(SQRT, makeMul(nomR) / makeMul(denR));
		}
	}
	else if (!denR.empty()) {
		den.emplace_back(SQRT, makeMul(denR));
	}

	if (f != 1.) nom.push_back(Symbolic(f));
	if (nom.empty()) nom.push_back(1.);
	ret = makeMul(nom);
	if (!ret.ahash()) return Symbolic(.0);
	if (!den.empty()) ret /= makeMul(den);

	return ret;
}

Symbolic pullConstant(const Symbolic& x_) {
	auto x = flattenMultiplications(x_);
	if (x.op() == MULC) return x;
	if (x.op() != MUL) return Symbolic(MULC, Symbolic(1.), x);

	double f = 1.;
	vector<Symbolic> childs;

	for (auto& c : x) {
		if (c.op() == CONST) f *= c.constant();
		else childs.push_back(c);
	}

	return Symbolic(MULC, Symbolic(f), makeSymbolic(MUL, childs));
}

Symbolic optimizeSum(const Symbolic& x0) {

	auto x = flattenAdditions(x0);
	if (x.op() != ADD) return x;

	map<Symbolic, double, AlgebraicHashFunctor> multiples;

	for (auto& c : mapData(x, pullConstant)) {
		assert(c.op() == MULC);
		multiples[c[1]] += c[0].constant();
	}

	vector<Symbolic> childs;
	for (auto& c : multiples) {
		if (c.second != .0 && c.second != -.0) {
			if (c.second == 1.) childs.push_back(c.first);
			else childs.push_back(flattenMultiplications(c.second * c.first));
		}
	}

	return makeSymbolic(ADD, childs);
}

/*
Symbolic reduceRoots(const Symbolic& x) {
	if(x.op() != SQRT) return x;
	assert(x.numChilds());

	auto r = x[0];
	if(r.op() == ADD) r = pullFactor(r);
	if(r.op() == MUL) r = flattenMultiplications(r);
	else return x;

	unordered_map<Symbolic, int, AlgebraicHashFunctor> factors;
	vector<Symbolic> out, in;

	for(auto& c : r) if(++factors[c] == 2) {
		out.push_back(c);
		factors[c] = 0;
	}

	for(auto& x : factors) {
		if(x.second == 1) in.push_back(x.first);
	}

	if(out.empty()) return x;
	auto ret = makeSymbolic(MUL, out);
	if(!in.empty()) {
		if(in.size() == 1 && in[0].op() == CONST)
			ret *= std::sqrt(in[0].constant());
		else ret *= Symbolic(SQRT, in);
	}

	return flattenMultiplications(ret);
}*/

// todo: flattenDivAndMul instead of flattenMultiplications
Symbolic pullFactor(Symbolic x) {

	if (x.op() != ADD) return x;
	if (x.numChilds() == 1) return x[0];
	if (x.numChilds() == 0) return .0;
	if (isSmallConstant(x.algebraicHash())) return Symbolic((double)x.algebraicHash());

	vector<vector<Symbolic>> summands;
	summands.reserve(x.numChilds());

	for (const auto& c : x) {
		if (c.op() == MUL || c.op() == MULC) {
			vector<Symbolic> factors;
			auto fc = flattenMultiplications(c);

			if (fc.op() == MUL) {// 'c' could be a multiplication with a single operand
				for (auto& cc : flattenMultiplications(c))
					factors.push_back(cc);
			}
			else factors.push_back(fc);

			summands.push_back(move(factors));

		}
		else summands.push_back({ c });
	}

	// find factor(s) that are shared by as many summands as possible
	map<Symbolic, set<int>, AlgebraicHashFunctor> factors;

	for (int i = 0; i < summands.size(); ++i) {
		for (auto& y : summands[i]) {
			factors[y].insert(i);
		}
	}

	// choose largest set of summands sharing a factor
	auto maxIt = max_element(factors.begin(), factors.end(), [&](const auto& a, const auto& b) {
		const auto scoreA = max((int)a.second.size() - 2, 0) * (a.first.complexity() + 1);
		const auto scoreB = max((int)b.second.size() - 2, 0) * (b.first.complexity() + 1);
		return scoreA < scoreB;
		});

	if (maxIt->second.size() < 2) return x; // no factor is present in more than one term

	auto it = maxIt->second.begin();
	vector<Symbolic> factored, rest;

	for (int i = 0; i < summands.size(); ++i) {
		if (it != maxIt->second.end() && i == *it) {

			auto its = find(summands[i].begin(), summands[i].end(), maxIt->first);
			assert(its != summands[i].end());
			summands[i].erase(its);

			if (summands[i].empty())
				factored.push_back(Symbolic(1.));
			else factored.emplace_back(MUL, summands[i]);

			++it;

		}
		else {
			rest.emplace_back(MUL, summands[i]);
		}
	}

	rest.push_back(maxIt->first * pullFactor(makeSymbolic(ADD, factored)));
	return pullFactor(flattenAdditions(makeSymbolic(ADD, rest)));
}

Symbolic findNegative(const Symbolic& expr) {
	if (expr.op() == MULC && expr[0].constant() == -1) return Symbolic(NEG, expr[1]);

	if (expr.op() == MUL) {
		vector<Symbolic> childs;
		int cnt = 0;

		for (auto& c : expr) {
			if (c.op() == CONST && c.constant() == -1.) ++cnt;
			else childs.push_back(c);
		}

		if (cnt % 2) return Symbolic(NEG, makeMul(childs));
		else return makeMul(childs);
	}
	else if (expr.op() == ADD) {

		vector<Symbolic> positive;
		vector<Symbolic> negative;

		for (auto& c : expr) {
			if (c.op() == NEG) negative.push_back(c[0]);
			else if (c.op() == MULC && c[0].constant() == -1.) negative.push_back(c[1]);
			else positive.push_back(c);
		}

		if (!negative.empty()) positive.emplace_back(NEG, makeAddition(negative));
		return makeAddition(positive);
	}

	return expr;
}

Symbolic simplifyExpression(const Symbolic& expr) {

	if (isSmallConstant(expr.ahash())) return (double)expr.ahash();
	if (expr.numChilds() == 0) return expr;

	auto x = optimizeSum(expr);
	x = pullFactor(x);
	x = optimizeProducts(x);
	return x;
}

Symbolic simplify(const Symbolic& expr) {

	return traverseGenerate<Symbolic>(expr,
		[](const Symbolic& x) {return simplifyExpression(x); },
		[](const Symbolic& x, const std::vector<Symbolic>& childs) {

			if (x.numChilds() == 0) return x;
			return findNegative(simplifyExpression(makeSymbolic(x.op(), childs)));

		}, true);

}

Symbolic removeConstantExpressions(const Symbolic& x) {
    return traverseGenerate<Symbolic>(x,
        [](const Symbolic& x) {
        return x;
            if(isSmallConstant(x.ahash())) {
                return Symbolic((double)x.ahash());
            } else return x;
        }, [](const Symbolic& x, const std::vector<Symbolic>& childs) {

            if (x.numChilds() == 0) return x;
            return makeSymbolic(x.op(), childs);

        }, true);
}

Symbolic referenceRedundant(Symbolic& x) {
	return traverseGenerate<Symbolic>(x,
		[](const Symbolic& x) {return x; },
		[](const Symbolic& x, const std::vector<Symbolic>& childs) {

			if (x.numChilds() == 0) return x;
			return makeSymbolic(x.op(), childs);

		}, true);
}

Symbolic referenceRedundant2(const Symbolic& x) {
	unordered_map<hash_t, Symbolic, IdentityHash<hash_t>> cache;
	vector<Symbolic> stack;

	prePostOrderTraverse(x, [&](const auto& y) {
		if (y.numChilds() == 0) {
			stack.push_back(y);
			return false;
		}

		if (isSmallConstant(y.ahash())) {
			stack.push_back(Symbolic((double)y.ahash()));
			return false;
		}

		auto it = cache.find(y.ahash());
		if (it != cache.end()) {
			stack.push_back(it->second);
			return false;
		}
		else return true;
		}, [&](const auto& y) {
			auto y2 = Symbolic(y.op(), std::vector<Symbolic>(stack.end() - y.numChilds(), stack.end()));
			stack.erase(stack.end() - y.numChilds(), stack.end());

			auto it = cache.find(y2.ahash());
			if (it != cache.end()) {
				y2 = it->second;
			}
			else {
				cache[y2.ahash()] = y2;
			}

			stack.push_back(y2);
		});

	assert(stack.size() == 1);
	return stack.front();
}

}
