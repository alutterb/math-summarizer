# Unconstrained Optimization 

## Contents

1.1. Optimality Conditions ..... p. 2
1.2. Gradient Methods - Convergence ..... p. 18
1.3. Gradient Methods - Rate of Convergence ..... p. 54
1.4. Newton's Method and Variations ..... p. 79
1.5. Least Squares Problems ..... p. 92
1.6. Conjugate Direction Methods ..... p. 118
1.7. Quasi-Newton Methods ..... p. 134
1.8. Nonderivative Methods ..... p. 142
1.9. Discrete-Time Optimal Control* ..... p. 149
1.10. Some Practical Guidelines ..... p. 165
1.11. Notes and Sources ..... p. 170

---

Mathematical models of optimization can be generally represented by a constraint set $X$ and a cost function $f$ that maps elements of $X$ into real numbers. The set $X$ consists of the available decisions $X$ and the cost $f(x)$ is a scalar measure of undesirability of choosing decision $x$. We want to find an optimal decision, that is, an $x^{*} \in X$ such that

$$
f\left(x^{*}\right) \leq f(x), \quad \forall x \in X
$$

In this book we focus on the case where each decision $x$ is an $n$-dimensional vector; that is, $x$ is an $n$-tuple of real numbers $\left(x_{1}, \ldots, x_{n}\right)$. Thus the constraint set $X$ is a subset of $\Re^{n}$, the $n$-dimensional Euclidean space. In this chapter we consider the unconstrained optimization problem where $X=\Re^{n}$, that is,

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } x \in \Re^{n}
\end{aligned}
$$

In subsequent chapters, we focus on problems where $X$ is specified by equality and inequality constraints.

For the most part, we assume in this chapter, that $f$ is a continuously differentiable function. In the first section, we discuss the main necessary and sufficient conditions for an optimal solution, and in the subsequent sections we consider algorithms for numerical computation of approximately optimal solutions.

# 1.1 OPTIMALITY CONDITIONS 

### 1.1.1 Variational Ideas

The main ideas underlying optimality conditions in nonlinear programming usually admit simple explanations although their detailed proofs are sometimes tedious. For this reason, we have chosen to first discuss informally these ideas in the present subsection, and to leave detailed statements of results and proofs for the next subsection.

## Local and Global Minima

A vector $x^{*}$ is an unconstrained local minimum of $f$ if it is no worse than its neighbors, that is, if there exists an $\epsilon>0$ such that

$$
f\left(x^{*}\right) \leq f(x), \quad \forall x \text { with }\left\|x-x^{*}\right\|<\epsilon
$$

---

(Unless stated otherwise, we use the standard Euclidean norm $\|x\|=\sqrt{x^{\prime} x}$. Appendix A describes in detail our mathematical notation and terminology.)

A vector $x^{*}$ is an unconstrained global minimum of $f$ if it is no worse than all other vectors, that is,

$$
f\left(x^{*}\right) \leq f(x), \quad \forall x \in \Re^{n}
$$

The unconstrained local or global minimum $x^{*}$ is said to be strict if the corresponding inequality above is strict for $x \neq x^{*}$. Figure 1.1.1 illustrates these definitions.
![[chapter_21_p3_img1.jpeg]]

Figure 1.1.1. Unconstrained local and global minima in one dimension.

The definitions of local and global minima can be extended to the case where $f$ is a function defined over a subset $C$ of $\Re^{n}$. In particular, we say that $x^{*}$ is a local minimum of $f$ over $C$ if $x^{*} \in C$ and there is an $\epsilon>0$ such that $f\left(x^{*}\right) \leq f(x)$ for all $x \in C$ with $\left\|x-x^{*}\right\|<\epsilon$. The definitions of a global and a strict minimum of $f$ over $C$ are analogous.

Local and global maxima are similarly defined. In particular, $x^{*}$ is an unconstrained local (global) maximum of $f$, if $x^{*}$ is an unconstrained local (global) minimum of the function $-f$.

# Necessary Conditions for Optimality 

If the cost function is differentiable, we can use gradients and [[Taylor Series|Taylor series expansions]] to compare the cost of a vector with the cost of its close neighbors. In particular, we consider small variations  from a given vector , which approximately, up to first order, yield a cost variation

$$
f\left(x^{*}+\Delta x\right)-f\left(x^{*}\right) \approx \nabla f\left(x^{*}\right)^{\prime} \Delta x
$$
[[Note 1 - Derivation of cost variation|Derivation HERE]]

---

and, up to second order, yield a cost variation

$$
f\left(x^{*}+\Delta x\right)-f\left(x^{*}\right) \approx \nabla f\left(x^{*}\right)^{\prime} \Delta x+\frac{1}{2} \Delta x^{\prime} \nabla^{2} f\left(x^{*}\right) \Delta x
$$

We expect that if $x^{*}$ is an unconstrained local minimum, the first order cost variation due to a small variation $\Delta x$ is nonnegative:

$$
\nabla f\left(x^{*}\right)^{\prime} \Delta x=\sum_{i=1}^{n} \frac{\partial f\left(x^{*}\right)}{\partial x_{i}} \Delta x_{i} \geq 0
$$

In particular, by taking $\Delta x$ to be [[Note 2 - Verbage Explanation|positive and negative multiples of the unit coordinate vectors]] (all coordinates equal to zero except for one which is equal to unity), we obtain $\partial f\left(x^{*}\right) / \partial x_{i} \geq 0$ and $\partial f\left(x^{*}\right) / \partial x_{i} \leq 0$, respectively, for all coordinates $i=1, \ldots, n$. Equivalently, we have the necessary condition (originally formulated by Fermat in 1629)

$$
\nabla f\left(x^{*}\right)=0
$$

This condition is proved in Prop. 1.1.1, given in the next subsection.
We also expect that the second order cost variation due to a small variation $\Delta x$ must also be nonnegative:

$$
\nabla f\left(x^{*}\right)^{\prime} \Delta x+\frac{1}{2} \Delta x^{\prime} \nabla^{2} f\left(x^{*}\right) \Delta x \geq 0
$$

Since $\nabla f\left(x^{*}\right)^{\prime} \Delta x=0$, we obtain

$$
\Delta x^{\prime} \nabla^{2} f\left(x^{*}\right) \Delta x \geq 0
$$

which implies that

$$
\nabla^{2} f\left(x^{*}\right): \text { positive semidefinite. }
$$

We prove this necessary condition in the subsequent Prop. 1.1.1. Appendix A reviews the definition and properties of positive definite and positive semidefinite matrices.

In what follows, we refer to a vector $x^{*}$ satisfying the condition $\nabla f\left(x^{*}\right)=0$ as a stationary point.

# The Case of a Convex Cost Function 

Convexity notions, reviewed in Appendix B, play a very important role in nonlinear programming. One reason is that when the cost function $f$ is convex, there is no distinction between local and global minima; every local minimum is also global. The idea is illustrated in Fig. 1.1.2 and the formal proof is given in Prop. B. 10 of Appendix B.

Another important fact is that the first order condition $\nabla f\left(x^{*}\right)=0$ is also sufficient for optimality if $f$ is convex. This is established in Prop. 1.1.2 in the next subsection. The proof is based on a basic property of a convex function $f$ : the linear approximation at a point $x^{*}$ based on the gradient, that is, $f\left(x^{*}\right)+\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right)$, underestimates $f(x)$, so if $\nabla f\left(x^{*}\right)=0$, then $f\left(x^{*}\right) \leq f(x)$ for all $x$ (see Prop. B. 3 in Appendix B).

Figure 1.1.3 shows how the first and second order necessary conditions can fail to guarantee local optimality of $x^{*}$ if $f$ is not convex.

---

![[chapter_21_p5_img2.jpeg]]

Figure 1.1.2. Illustration of why local minima of convex functions are also global. Suppose that $f$ is convex and that $x^{*}$ is a local minimum of $f$. Let $\bar{x}$ be such that $f(\bar{x})<f\left(x^{*}\right)$. By convexity, for all $\alpha \in(0,1)$,

$$
f\left(\alpha x^{*}+(1-\alpha) \bar{x}\right) \leq \alpha f\left(x^{*}\right)+(1-\alpha) f(\bar{x})<f\left(x^{*}\right)
$$

Thus, $f$ decreases monotonically on the line segment connecting $x^{*}$ with $\bar{x}$ and $x^{*}$ cannot be a local minimum which is not global.
![[chapter_21_p5_img3.jpeg]]

Figure 1.1.3. Illustration of the first order necessary optimality condition of zero slope $\left[\nabla f\left(x^{*}\right)=0\right]$ and the second order necessary optimality condition of nonnegative curvature $\left[\nabla^{2} f\left(x^{*}\right) \geq 0\right]$ for functions of one variable. The first order condition is satisfied not only by local minima, but also by local maxima and "inflection" points, such as the one on the middle figure above. In some cases [e.g. the functions $f(x)=x^{3}$ and $f(x)=-|x|^{3}$ ] the second order condition is also satisfied by local maxima and inflection points. If the function $f$ is convex, the condition $\nabla f\left(x^{*}\right)=0$ is necessary and sufficient for global optimality of $x^{*}$.

---

# Sufficient Conditions for Optimality 

Suppose we have a vector $x^{*}$ that satisfies the first order necessary optimality condition

$$
\nabla f\left(x^{*}\right)=0
$$

and also satisfies the following strengthened form of the second order necessary optimality condition

$$
\nabla^{2} f\left(x^{*}\right): \text { positive definite }
$$

(that is, the Hessian is positive definite rather than semidefinite). Then, for all $\Delta x \neq 0$ we have

$$
\Delta x^{\prime} \nabla^{2} f\left(x^{*}\right) \Delta x>0
$$

implying that at $x^{*}$ the second order variation of $f$ due to a small nonzero variation $\Delta x$ is positive. Thus, $f$ tends to increase strictly with small excursions from $x^{*}$, suggesting that the above conditions (1.1) and (1.2) are sufficient for local optimality of $x^{*}$. This is established in Prop. 1.1.3.

Local minima that don't satisfy the sufficiency conditions (1.1) and (1.2) are called singular; otherwise they are called nonsingular. Singular local minima are harder to deal with for two reasons. First, in the absence of convexity of $f$, their optimality cannot be ascertained using easily verifiable sufficiency conditions. Second, in their neighborhood, the behavior of the most commonly used optimization algorithms tends to be slow and/or erratic, as we will see in the subsequent sections.

## Quadratic Cost Functions

Consider the quadratic function

$$
f(x)=\frac{1}{2} x^{\prime} Q x-b^{\prime} x
$$

where $Q$ is a symmetric $n \times n$ matrix and $b$ is a vector in $\Re^{n}$. If $x^{*}$ is a local minimum of $f$, we must have, by the necessary optimality conditions,

$$
\nabla f\left(x^{*}\right)=Q x^{*}-b=0, \quad \nabla^{2} f\left(x^{*}\right)=Q: \text { positive semidefinite. }
$$

Thus, if $Q$ is not positive semidefinite, $f$ can have no local minima. If $Q$ is positive semidefinite, $f$ is convex [Prop. B.4(d) of Appendix B], so any vector $x^{*}$ satisfying the first order condition $\nabla f\left(x^{*}\right)=Q x^{*}-b=0$ is a global minimum of $f$. On the other hand there might not exist a solution of the equation $\nabla f\left(x^{*}\right)=Q x^{*}-b=0$ if $Q$ is singular. If, however, $Q$ is positive definite (and hence invertible, by Prop. A. 20 of Appendix A), the equation $Q x^{*}-b=0$ can be solved uniquely and the vector $x^{*}=Q^{-1} b$

---

is the unique global minimum. This is consistent with Prop. 1.1.2(a) to be given shortly, which asserts that strictly convex functions can have at most one global minimum $\mid f$ is strictly convex if and only if $Q$ is positive definite; Prop. B.4(e) of Appendix B]. Figure 1.1.4 illustrates the various special cases considered.

Quadratic cost functions are important in nonlinear programming because they arise frequently in applications, but they are also important for another reason. From the Taylor expansion

$$
f(x)=f\left(x^{*}\right)+\frac{1}{2}\left(x-x^{*}\right)^{\prime} \nabla^{2} f\left(x^{*}\right)\left(x-x^{*}\right)+o\left(\left\|x-x^{*}\right\|^{2}\right)
$$

it is seen that a nonquadratic cost function can be approximated well by a quadratic function near a nonsingular local minimum $x^{*} \mid \nabla^{2} f\left(x^{*}\right)$ : positive definite]. This means that we can carry out much of our analysis and experimentation with algorithms using positive definite quadratic functions and expect that the conclusions will largely carry over to more general cost functions near convergence to such local minima. However, for local minima near which the Hessian matrix either does not exist or is singular, the higher than second order terms in the Taylor expansion are not negligible and an algorithmic analysis based on quadratic cost functions will likely be seriously flawed.

# Existence of Optimal Solutions 

In many cases it is important to know that there exists at least one global minimum of a function $f$ over a set $X$. Generally, such a minimum need not exist. For example, the scalar functions $f(x)=x$ and $f(x)=e^{x}$ have no global minima over the set of real numbers. The first function decreases without bound to $-\infty$ as $x$ tends toward $-\infty$, while the second decreases toward 0 as $x$ tends toward $-\infty$ but always takes positive values. Given the range of values that $f(x)$ takes as $x$ ranges over $X$, that is, the set of real numbers

$$
\{f(x) \mid x \in X\}
$$

there are two possibilities:

1. The set $\{f(x) \mid x \in X\}$ is bounded below; that is, there exists a scalar $M$ such that $M \leq f(x)$ for all $x \in X$. In this case, the greatest lower bound of $\{f(x) \mid x \in X\}$ is a real number, which is denoted by $\inf _{x \in X} f(x)$. For example, $\inf _{x \in \Re} e^{x}=0$ and $\inf _{x<0} e^{x}=0$.
2. The set $\{f(x) \mid x \in S\}$ is unbounded below (i.e., contains arbitrarily small real numbers). In this case we write

$$
\inf _{x \in X} f(x)=-\infty
$$

Existence of at least one global minimum is guaranteed if $f$ is a continuous function and $X$ is a compact subset of $\Re^{n}$. This is the Weierstrass

---

![[chapter_21_p8_img4.jpeg]]

Figure 1.1.4. Illustration of the isocost surfaces of the quadratic cost function $f: \Re^{2} \mapsto \Re$ given by

$$
f(x, y)=\frac{1}{2}\left(\alpha x^{2}+\beta y^{2}\right)-x
$$

for various values of $\alpha$ and $\beta$.
theorem. (see Prop. A. 8 in Appendix A). By a related result, also shown in Prop. A. 8 of Appendix A, existence of an optimal solution is guaranteed if $f: \Re^{n} \mapsto \Re$ is a continuous function, $X$ is closed, and $f$ is coercive, that is, $f(x) \rightarrow \infty$ when $\|x\| \rightarrow \infty$.

# Why do we Need Optimality Conditions? 

Hardly anyone would doubt that optimality conditions are fundamental to the analysis of an optimization problem. Yet in practice optimality

---

conditions are used in different ways than one might speculate.

The most straightforward method to use optimality conditions to solve an optimization problem, is as follows: First, find all points satisfying the first order necessary condition $\nabla f(x)=0$; then (if $f$ is not known to be convex), check the second order necessary condition ( $\nabla^{2} f$ : positive semidefinite) for each of these points, filtering out those that do not satisfy it; finally for the remaining candidates, check if $\nabla^{2} f$ is positive definite, in which case we are sure that they are strict local minima.

A slightly different alternative method is to find all points satisfying the necessary conditions, and to declare as global minimum the one with smallest cost value. However, here it is essential to know that a global minimum exists. As an example, for the one-dimensional function $f(x)=$ $x^{2}-x^{4}$, the points satisfying the necessary condition $\nabla f(x)=2 x-4 x^{3}=0$ are $0,1 / \sqrt{2}$, and $-1 / \sqrt{2}$, and of these, 0 gives the smallest cost value. However, we cannot declare 0 as the global minimum, because we don't know if a global minimum exists. Indeed, in this example none of the points $0,1 / \sqrt{2}$, and $-1 / \sqrt{2}$ is a global minimum, because $f$ decreases to $-\infty$ as $|x| \rightarrow \infty$, and has no global minimum.

It is important to realize, however, that except under very favorable circumstances (usually for one- or two-dimensional problems), using optimality conditions as described above does not work. The reason is that solving for $x$ the system of equations $\nabla f(x)=0$ is nontrivial; algorithmically, it is usually as difficult as solving the original optimization problem.

The principal context in which optimality conditions become useful will not become apparent until we consider iterative optimization algorithms in subsequent sections. We will see that optimality conditions often provide the basis for the development and the analysis of algorithms. In particular, algorithms recognize solutions by checking whether they satisfy various optimality conditions and terminate when such conditions hold approximately. Furthermore, the behavior of various algorithms in the neighborhood of a local minimum often depends on whether various optimality conditions are satisfied at that minimum. Thus, for example, sufficiency conditions play a key role in assertions regarding the speed of convergence of various algorithms.

There is one other important context, prominently arising in microeconomic theory, where optimality conditions provide the basis for analysis. Here one is interested primarily not in finding an optimal solution, but rather in how the optimal solution is affected by changes in the problem data. For example, an economist may be interested in how the prices of some raw materials will affect the availability of certain goods that are produced by using these raw materials; the assumption here is that the amounts produced are the variables of a profit optimization problem, which is solved by the corresponding producers. This type of analysis is known as sensitivity analysis, and is discussed next.

---

# Sensitivity* 

Suppose that we want to quantify the variation of the optimal solution as a vector of parameters changes. In particular, consider the optimization problem

$$
\begin{aligned}
& \operatorname{minimize} f(x, a) \\
& \text { subject to } x \in \Re^{n}
\end{aligned}
$$

where $f: \Re^{m+n} \mapsto \Re$ is a twice continuously differentiable function involving the $m$-dimensional parameter vector $a$. Let $x(a)$ denote the global minimum corresponding to $a$, assuming for the moment that it exists and is unique. By the first order necessary condition we have

$$
\nabla_{x} f(x(a), a)=0, \quad \forall a \in \Re^{m}
$$

and by differentiating this relation with respect to $a$, we obtain

$$
\nabla x(a) \nabla_{x x}^{2} f(x(a), a)+\nabla_{x a} f(x(a), a)=0
$$

where the elements of the $m \times n$ gradient matrix $\nabla x(a)$ are the first partial derivatives of the components of $x(a)$ with respect to the different components of $a$. Assuming that the inverse below exists, we have

$$
\nabla x(a)=-\nabla_{x a} f(x(a), a)\left(\nabla_{x x}^{2} f(x(a), a)\right)^{-1}
$$

which gives the first order variation of the components of the optimal $x$ with respect to the components of $a$.

For the preceding analysis to be precise, we must be sure that $x(a)$ exists and is differentiable as a function of $a$. The principal analytical framework for this is the implicit function theorem (Prop. A. 25 in Appendix A). With the aid of this theorem, we can define $x(a)$ in some sphere around a minimum $\bar{x}=x(\bar{a})$ corresponding to a nominal parameter value $\bar{a}$, assuming that the Hessian matrix $\nabla_{x x}^{2} f(\bar{x}, \bar{a})$ is positive definite. Thus, the preceding development and the formula (1.5) for the matrix $\nabla x(a)$ can be justified provided the nominal local minimum $\bar{x}$ is nonsingular.

We will postpone further discussion of sensitivity analysis for Section 3.2 , where we will provide a constrained version of the expression (1.5) for $\nabla x(a)$.

### 1.1.2 Main Optimality Conditions

We now provide formal statements and proofs of the optimality conditions.

---

Proposition 1.1.1: (Necessary Optimality Conditions) Let $x^{*}$ be an unconstrained local minimum of $f: \Re^{n} \mapsto \Re$, and assume that $f$ is continuously differentiable in an open set $S$ containing $x^{*}$. Then

$$
\nabla f\left(x^{*}\right)=0 . \quad \text { (First Order Necessary Condition) }
$$

If in addition $f$ is twice continuously differentiable within $S$, then
$\nabla^{2} f\left(x^{*}\right)$ : positive semidefinite. (Second Order Necessary Condition)

Proof: Fix some $d \in \Re^{n}$. Then, using the chain rule to differentiate the function $g(\alpha)=f\left(x^{*}+\alpha d\right)$ of the scalar $\alpha$, we have

$$
0 \leq \lim _{\alpha \downarrow 0} \frac{f\left(x^{*}+\alpha d\right)-f\left(x^{*}\right)}{\alpha}=\frac{d g(0)}{d \alpha}=d^{\prime} \nabla f\left(x^{*}\right)
$$

where the inequality follows from the assumption that $x^{*}$ is a local minimum. Since $d$ is arbitrary, the same inequality holds with $d$ replaced by $-d$. Therefore, $d^{\prime} \nabla f\left(x^{*}\right)=0$ for all $d \in \Re^{n}$, which shows that $\nabla f\left(x^{*}\right)=0$.

Assume that $f$ is twice continuously differentiable, and let $d$ be any vector in $\Re^{n}$. For all $\alpha \in \Re$, the second order Taylor series expansion yields

$$
f\left(x^{*}+\alpha d\right)-f\left(x^{*}\right)=\alpha \nabla f\left(x^{*}\right)^{\prime} d+\frac{\alpha^{2}}{2} d^{\prime} \nabla^{2} f\left(x^{*}\right) d+o\left(\alpha^{2}\right)
$$

Using the condition $\nabla f\left(x^{*}\right)=0$ and the local optimality of $x^{*}$, we see that there is a sufficiently small $\epsilon>0$ such that for all $\alpha$ with $\alpha \in(0, \epsilon)$,

$$
0 \leq \frac{f\left(x^{*}+\alpha d\right)-f\left(x^{*}\right)}{\alpha^{2}}=\frac{1}{2} d^{\prime} \nabla^{2} f\left(x^{*}\right) d+\frac{o\left(\alpha^{2}\right)}{\alpha^{2}}
$$

Taking the limit as $\alpha \rightarrow 0$ and using the fact

$$
\lim _{\alpha \rightarrow 0} \frac{o\left(\alpha^{2}\right)}{\alpha^{2}}=0
$$

we obtain $d^{\prime} \nabla^{2} f\left(x^{*}\right) d \geq 0$, showing that $\nabla^{2} f\left(x^{*}\right)$ is positive semidefinite. Q.E.D.

The following proposition handles the case of a convex cost.

---

Proposition 1.1.2: (Convex Cost Function) Let $f: C \mapsto \Re$ be a convex function over the convex set $C$.
(a) A local minimum of $f$ over $C$ is also a global minimum over $C$. If in addition $f$ is strictly convex, then there exists at most one global minimum of $f$.
(b) If $f$ is convex and the set $C$ is open, then $\nabla f\left(x^{*}\right)=0$ is a necessary and sufficient condition for a vector $x^{*} \in C$ to be a global minimum of $f$ over $C$.

Proof: Part (a) is proved in Prop. B. 10 of Appendix B. To show part (b), note that by Prop. B. 3 of Appendix B, we have

$$
f(x) \geq f\left(x^{*}\right)+\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right), \quad \forall x \in C
$$

If $\nabla f\left(x^{*}\right)=0$, we obtain $f(x) \geq f\left(x^{*}\right)$ for all $x \in C$, so $x^{*}$ is a global minimum. Q.E.D.

In the absence of convexity, we have the following sufficiency conditions for local optimality.

Proposition 1.1.3: (Second Order Sufficient Optimality Conditions) Let $f: \Re^{n} \mapsto \Re$ be twice continuously differentiable in an open set $S$. Suppose that a vector $x^{*} \in S$ satisfies the conditions

$$
\nabla f\left(x^{*}\right)=0, \quad \nabla^{2} f\left(x^{*}\right): \text { positive definite. }
$$

Then, $x^{*}$ is a strict unconstrained local minimum of $f$. In particular, there exist scalars $\gamma>0$ and $\epsilon>0$ such that

$$
f(x) \geq f\left(x^{*}\right)+\frac{\gamma}{2}\left\|x-x^{*}\right\|^{2}, \quad \forall x \text { with }\left\|x-x^{*}\right\|<\epsilon
$$

Proof: Let $\lambda$ be the smallest eigenvalue of $\nabla^{2} f\left(x^{*}\right)$. By Prop. A. $20(\mathrm{~b})$ of Appendix A, $\lambda$ is positive since $\nabla^{2} f\left(x^{*}\right)$ is positive definite. Furthermore, by Prop. A. $18(\mathrm{~b})$ of Appendix A, $d^{\prime} \nabla^{2} f\left(x^{*}\right) d \geq \lambda\|d\|^{2}$ for all $d \in \Re^{n}$. Using this relation, the hypothesis $\nabla f\left(x^{*}\right)=0$, and the second order Taylor series expansion, we have for all $d$

$$
\begin{aligned}
f\left(x^{*}+d\right)-f\left(x^{*}\right) & =\nabla f\left(x^{*}\right)^{\prime} d+\frac{1}{2} d^{\prime} \nabla^{2} f\left(x^{*}\right) d+o\left(\|d\|^{2}\right) \\
& \geq \frac{\lambda}{2}\|d\|^{2}+o\left(\|d\|^{2}\right)
\end{aligned}
$$

---

$$
=\left(\frac{\lambda}{2}+\frac{o\left(\\|d\|^{2}\right)}{\\|d\|^{2}}\right)\|d\|^{2}
$$

It is seen that Eq. (1.6) is satisfied for any $\epsilon>0$ and $\gamma>0$ such that

$$
\frac{\lambda}{2}+\frac{o\left(\\|d\|^{2}\right)}{\\|d\|^{2}} \geq \gamma, \quad \forall d \text { with }\|d\|<\epsilon
$$

Q.E.D.

# E XERCISES 

## 1.1

For each value of the scalar $\beta$, find the set of all stationary points of the following function of the two variables $x$ and $y$

$$
f(x, y)=x^{2}+y^{2}+\beta x y+x+2 y
$$

Which of these stationary points are local minima and which are global minima?

## 1.2

In each of the following problems fully justify your answer using optimality conditions.
(a) Show that the 2-dimensional function $f(x, y)=\left(x^{2}-4\right)^{2}+y^{2}$ has two global minima and one stationary point, which is neither a local maximum nor a local minimum.
(b) Show that the 2-dimensional function $f(x, y)=\left(y-x^{2}\right)^{2}-x^{2}$ has only one stationary point, which is neither a local maximum nor a local minimum.
(c) Find all local minima of the 2-dimensional function $f(x, y)=\frac{1}{2} x^{2}+$ $x \cos y$.
(d) Find all local minima and all local maxima of the 2-dimensional function $f(x, y)=\sin x+\sin y+\sin (x+y)$ within the set $\{(x, y) \mid 0<x<2 \pi, 0<$ $y<2 \pi\}$.

---

# 1.3 [Hes75] 

Let $f: \Re^{n} \mapsto \Re$ be a differentiable function. Suppose that a point $x^{*}$ is a local minimum of $f$ along every line that passes through $x^{*}$; that is, the function

$$
g(\alpha)=f\left(x^{*}+\alpha d\right)
$$

is minimized at $\alpha=0$ for all $d \in \Re^{n}$.
(a) Show that $\nabla f\left(x^{*}\right)=0$.
(b) Show by example that $x^{*}$ need not be a local minimum of $f$. Hint: Consider the function of two variables $f(y, z)=\left(z-p y^{2}\right)\left(z-q y^{2}\right)$, where $0<p<q$; see Fig. 1.1.5. Show that $(0,0)$ is a local minimum of $f$ along every line that passes through $(0,0)$. Furthermore, if $p<m<q$, then $f\left(y, m y^{2}\right)<0$ if $y \neq 0$ while $f(0,0)=0$.
![[chapter_21_p14_img5.jpeg]]

Figure 1.1.5. Three-dimensional graph of the function $f(y, z)=\left(z-p y^{2}\right)\left(z-q y^{2}\right)$ for $p=1$ and $q=4$ (cf. Exercise 1.3). The origin is a local minimum with respect to every line that passes through it, but is not a local minimum of $f$.

## 1.4

Use optimality conditions to show that for all $x>0$ we have

$$
\frac{1}{x}+x \geq 2
$$

---

1.5

Find the parallelepiped of unit volume that has the minimum surface area. Hint: By eliminating one of the dimensions, show that the problem is equivalent to the minimization over $x>0$ and $y>0$ of

$$
f(x, y)=x y+\frac{1}{x}+\frac{1}{y}
$$

# 1.6 (The Weber Point of a Set of Points) 

We want to find a point $x$ in the plane whose sum of weighted distances from a given set of points $y_{1}, \ldots, y_{m}$ is minimized. Mathematically, the problem is

$$
\begin{aligned}
& \operatorname{minimize} \sum_{i=1}^{m} w_{i}\left\|x-y_{i}\right\| \\
& \text { subject to } x \in \Re^{n}
\end{aligned}
$$

where $w_{1}, \ldots, w_{m}$ are given positive scalars. Show that there exists a global minimum for this problem and that it can be realized by means of the mechanical model shown in Fig. 1.1.6. Is the optimal solution always unique?
![[chapter_21_p15_img6.jpeg]]

Figure 1.1.6. Mechanical model (known as the Varignon frame) associated with the Weber problem (Exercise 1.6). It consists of a board with a hole drilled at each of the given points $y_{i}$. Through each hole, a string is passed with the corresponding weight $w_{i}$ attached. The other ends of the strings are tied with a knot as shown. In the absence of friction or tangled strings, the forces at the knot reach equilibrium when the knot is located at an optimal solution $x^{*}$.

---

# 1.7 (Diffraction Law) 

Let $p$ and $q$ be two points on the plane that lie on opposite sides of a horizontal axis. Assume that the speed of light from $p$ and from $q$ to the horizontal axis is $v$ and $w$, respectively, and that light reaches a point from other points along paths of minimum travel time. Find the path that a ray of light would follow from $p$ to $q$.

## 1.8 (Steiner's Problem)

Given a triangle in the plane, consider the problem of finding a point whose sum of distances from the vertices of the triangle is minimal. Show that such a point is either a vertex, or else it is such that each side of the triangle is seen from that point at a 120 degree angle (this is known as the Torricelli point).

## 1.9 (Stability)

We are often interested in whether optimal solutions change radically when the problem data are slightly perturbed. This issue is addressed by stability analysis, to be contrasted with sensitivity analysis, which deals with how much optimal solutions change when problem data change. An unconstrained local minimum $x^{*}$ of a function $f$ is said to be locally stable if there exists a $\delta>0$ such that all sequences $\left\{x^{k}\right\}$ with $f\left(x^{k}\right) \rightarrow f\left(x^{*}\right)$ and $\left\|x^{k}-x^{*}\right\|<\delta$, for all $k$, converge to $x^{*}$. (For example, a strict local minimum of a continuous function is locally stable.) Suppose that $x^{*}$ is a locally stable local minimum of the continuous function $f$ and let $g$ be a continuous function. Show that there exists a $\delta>0$ such that for all sufficiently small $\epsilon>0$, the function $f(x)+\epsilon g(x)$ has an unconstrained local minimum $x_{\epsilon}$ that lies within the sphere centered at $x^{*}$ with radius $\delta$. Furthermore, $x_{\epsilon} \rightarrow x^{*}$ as $\epsilon \rightarrow 0$.

### 1.10 (Nonconvex Level sets [Dun87])

Let $f: \Re^{2} \mapsto \Re$ be the function

$$
f(x)=x_{2}^{2}-a x_{2}\|x\|^{2}+\|x\|^{4}
$$

where $0<a<2$ (see Fig. 1.1.7). Show that $f(x)>0$ for all $x \neq 0$, so that the origin is the unique global minimum. Show also that there exists a $\bar{\gamma}>0$ such that for all $\gamma \in(0, \bar{\gamma}]$, the level set $L_{\gamma}=\{x \mid f(x) \leq \gamma\}$ is not convex. Hint: Show that for $\gamma \in(0, \bar{\gamma}]$, there is a $p>0$ and a $q>0$ such that the vectors $(-p, q)$ and $(p, q)$ belong to $L_{\gamma}$, but $(0, q)$ does not belong to $L_{\gamma}$.

---

![[chapter_21_p17_img7.jpeg]]

Figure 1.1.7. Level sets of the function $f$ of Exercise 1.10 for the case where $a=1.98$. The unique global minimum is the origin, but the level sets of $f$ are nonconvex.

# 1.11 (Singular Strict Local Minima [Dun87]) 

Show that if $x^{*}$ is a [[Vocab List#Nonsingular (nondegenerate) stationary point|nonsingular]] strict local minimum of $\alpha$ twice continuously differentiable function stationary, then $x^{*}$ is an [[Vocab List#Isolated stationary point|isolated stationary]] point; that is, there is a sphere centered at  such that $x^{*}$ is the only stationary point of  within that sphere. Use the following example function $f: \Re \mapsto \Re$ to show that this need not be true if $x^{*}$ is a [[Vocab List#Singular (degenerate) stationary point|singular]] strict local minimum:

$$
f(x)= \begin{cases}x^{2}\left(\sqrt{2}-\sin \left(\frac{\pi}{6}-\sqrt{3} \ln x^{2}\right)\right) & \text { if } x \neq 0 \\ 0 & \text { if } x=0\end{cases}
$$

In particular, show that $x^{*}=0$ is the unique (singular) global minimum, while the sequence $\left\{x^{k}\right\}$ of nonsingular local minima, where

$$
x^{k}=e^{\frac{(1-9 k) \pi}{\pi \sqrt{3}}}
$$

converges to $x^{*}$.

---

# 1.2 GRADIENT METHODS - CONVERGENCE 

We now start our development of computational methods for unconstrained optimization. The conceptual framework of this section is fundamental in nonlinear programming and applies to constrained optimization methods as well.

### 1.2.1 Descent Directions and Stepsize Rules

As in the case of optimality conditions, the main ideas of unconstrained optimization methods have simple geometrical explanations, but the corresponding convergence analysis is often complex. For this reason we first discuss informally the methods and their behavior in the present subsection, and we substantiate our conclusions with rigorous analysis in the next subsection.

Consider the problem of unconstrained minimization of a continuously differentiable function $f: \Re^{n} \mapsto \Re$. Most of the interesting algorithms for this problem rely on an important idea, called iterative descent that works as follows: We start at some point $x^{0}$ (an initial guess) and successively generate vectors $x^{1}, x^{2}, \ldots$, such that $f$ is decreased at each iteration, that is

$$
f\left(x^{k+1}\right)<f\left(x^{k}\right), \quad k=0,1, \ldots
$$

(cf. Fig. 1.2.1). In doing so, we successively improve our current solution estimate and we hope to decrease $f$ all the way to its minimum. In this section, we introduce a general class algorithms based on iterative descent, and we analyze their convergence to local minima. In the next section we examine their rate of convergence properties.
![[chapter_21_p18_img8.jpeg]]

Figure 1.2.1. Iterative descent for minimizing a function $f$. Each vector in the generated sequence has a lower cost than its predecessor.

---

# Gradient Methods 

Given a vector $x \in \Re^{n}$ with $\nabla f(x) \neq 0$, consider the half line of vectors

$$
x_{\alpha}=x-\alpha \nabla f(x), \quad \forall \alpha \geq 0
$$

From the first order Taylor series expansion around $x$ we have

$$
\begin{aligned}
f\left(x_{\alpha}\right) & =f(x)+\nabla f(x)^{\prime}\left(x_{\alpha}-x\right)+o\left(\left\|x_{\alpha}-x\right\|\right) \\
& =f(x)-\alpha\|\nabla f(x)\|^{2}+o(\alpha\|\nabla f(x)\|)
\end{aligned}
$$

so we can write

$$
f\left(x_{\alpha}\right)=f(x)-\alpha\|\nabla f(x)\|^{2}+o(\alpha)
$$

The term $\alpha\|\nabla f(x)\|^{2}$ dominates $o(\alpha)$ for $\alpha$ near zero, so for positive but sufficiently small $\alpha, f\left(x_{\alpha}\right)$ is smaller than $f(x)$ as illustrated in Fig. 1.2.2.
![[chapter_21_p19_img9.jpeg]]

Figure 1.2.2. If $\nabla f(x) \neq 0$, there is an interval $(0, \delta)$ of stepsizes such that $f(x-$ $\alpha \nabla f(x))<f(x)$ for all $\alpha \in(0, \delta)$.

Carrying this idea one step further, consider the half line of vectors

$$
x_{\alpha}=x+\alpha d, \quad \forall \alpha \geq 0
$$

where the direction vector $d \in \Re^{n}$ makes an angle with $\nabla f(x)$ that is greater than 90 degrees, that is,

$$
\nabla f(x)^{\prime} d<0
$$

Again by Taylor's Theorem we have

$$
f\left(x_{\alpha}\right)=f(x)+\alpha \nabla f(x)^{\prime} d+o(\alpha)
$$

---

![[chapter_21_p20_img10.jpeg]]

Figure 1.2.3. If the direction $d$ makes an angle with $\nabla f(x)$ that is greater than 90 degrees, that is, $\nabla f(x)^{\prime} d<0$, there is an interval $(0, \delta)$ of stepsizes such that $f(x+$ $\alpha d)<f(x)$ for all $\alpha \in(0, \delta)$.

For $\alpha$ near zero, the term $\alpha \nabla f(x)^{\prime} d$ dominates $o(\alpha)$ and as a result, for positive but sufficiently small $\alpha . f(x+\alpha d)$ is smaller than $f(x)$ as illustrated in Fig. 1.2.3.

The preceding observations form the basis for the broad and important class of algorithms

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}, \quad k=0,1, \ldots
$$

where, if $\nabla f\left(x^{k}\right) \neq 0$, the direction $d^{k}$ is chosen so that

$$
\nabla f\left(x^{k}\right)^{\prime} d^{k}<0
$$

and the stepsize $\alpha^{k}$ is chosen to be positive. If $\nabla f\left(x^{k}\right)=0$, then the method stops, that is $x^{k+1}=x^{k}$ (equivalently we choose $d^{k}=0$ ). In view of the relation (2.2) of the direction $d^{k}$ and the gradient $\nabla f\left(x^{k}\right)$, we call algorithms of this type gradient methods. [There is no universally accepted name for these algorithms; some authors reserve the name "gradient method" for the special case where $d^{k}=-\nabla f\left(x^{k}\right)$.] The great majority of the gradient methods that we will consider are also descent algorithms; that is, the stepsize $\alpha^{k}$ is selected so that

$$
f\left(x^{k}+\alpha^{k} d^{k}\right)<f\left(x^{k}\right), \quad k=0,1, \ldots
$$

However, there are some exceptions.
There is a large variety of possibilities for choosing the direction $d^{k}$ and the stepsize $\alpha^{k}$ in a gradient method. Indeed there is no single gradient method that can be recommended for all or even most problems. Otherwise said, given any one of the numerous methods and variations thereof that we will discuss, there is a class of problems for which this method is wellsuited. Our principal analytical aim is to develop a few guiding principles for understanding the performance of broad classes of methods and for appreciating the practical contexts in which their use is most appropriate.

---

# Selecting the Descent Direction 

Many gradient methods are specified in the form

$$
x^{k+1}=x^{k}-\alpha^{k} D^{k} \nabla f\left(x^{k}\right)
$$

where $D^{k}$ is a positive definite symmetric matrix. In this case, $d^{k}=$ $-D^{k} \nabla f\left(x^{k}\right)$, and the descent condition $\nabla f\left(x^{k}\right)^{\prime} d^{k}<0$, which can be written as

$$
\nabla f\left(x^{k}\right)^{\prime} D^{k} \nabla f\left(x^{k}\right)>0
$$

holds thanks to the positive definiteness of $D^{k}$.
Here are some examples of choices of the matrix $D^{k}$, resulting in methods that are widely used:

## Steepest Descent

$$
D^{k}=I . \quad k=0,1 \ldots
$$

where $I$ is the $n \times n$ identity matrix. This is the simplest choice but it often leads to slow convergence, as we will see in Section 1.3. The difficulty is illustrated in Fig. 1.2.4 and motivates the methods of the subsequent examples.
![[chapter_21_p21_img11.jpeg]]

Figure 1.2.4. Slow convergence of the steepest descent method

$$
x^{k+1}=x^{k}-\alpha^{k} \nabla f\left(x^{k}\right)
$$

when the equal cost surfaces of $f$ are "elongated." The difficulty is that the gradient direction is almost orthogonal to the direction that leads to the minimum. As a result the method is zig-zagging without making fast progress.

## Newton's Method

$$
D^{k}=\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1}, \quad k=0,1, \ldots
$$

provided $\nabla^{2} f\left(x^{k}\right)$ is positive definite. If $\nabla^{2} f\left(x^{k}\right)$ is not positive definite, some modification is necessary as will be explained in Section 1.4. The idea in Newton's method is to minimize at each iteration the quadratic approximation of $f$ around the current point $x^{k}$ given by

$$
f^{k}(x)=f\left(x^{k}\right)+\nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)+\frac{1}{2}\left(x-x^{k}\right)^{\prime} \nabla^{2} f\left(x^{k}\right)\left(x-x^{k}\right)
$$

---

(see Fig. 1.2.5). By setting the derivative of $f^{k}(x)$ to zero,

$$
\nabla f\left(x^{k}\right)+\nabla^{2} f\left(x^{k}\right)\left(x-x^{k}\right)=0
$$

we obtain the next iterate $x^{k+1}$ as the minimum of $f^{k}(x)$ :

$$
x^{k+1}=x^{k}-\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1} \nabla f\left(x^{k}\right)
$$

This is the pure Newton iteration. It corresponds to the more general iteration

$$
x^{k+1}=x^{k}-\alpha^{k}\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1} \nabla f\left(x^{k}\right)
$$

where the stepsize $\alpha^{k}=1$. Note that Newton's method finds the global minimum of a positive definite quadratic function in a single iteration (assuming $\alpha^{k}=1$ ). More generally, Newton's method typically converges very fast asymptotically and does not exhibit the zig-zagging behavior of steepest descent, as we will show in Section 1.4. For this reason many other methods try to emulate Newton's method. Some examples are given below.
![[chapter_21_p22_img12.jpeg]]

Figure 1.2.5. Illustration of the fast convergence rate of Newton's method with a stepsize equal to one. Given $x^{k}$, the method obtains $x^{k+1}$ as the minimum of a quadratic approximation of $f$ based on a second order Taylor expansion around $x^{k}$.

Diagonally Scaled Steepest Descent

$$
D^{k}=\left(\begin{array}{cccccccc}
d_{1}^{k} & 0 & 0 & \cdots & 0 & 0 & 0 \\
0 & d_{2}^{k} & 0 & \cdots & 0 & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \\
0 & 0 & 0 & \cdots & 0 & d_{n-1}^{k} & 0 \\
0 & 0 & 0 & \cdots & 0 & 0 & d_{n}^{k}
\end{array}\right) . \quad k=0,1, \ldots
$$

---

where $d_{i}^{k}$ are positive scalars, thus ensuring that $D^{k}$ is positive definite. A popular choice, resulting in a method known as a diagonal approximation to Newton's method, is to take $d_{i}^{k}$ to be an approximation to the inverted second partial derivative of $f$ with respect to $x_{i}$, that is,

$$
d_{i}^{k} \approx\left(\frac{\partial^{2} f\left(x^{k}\right)}{\left(\partial x_{i}\right)^{2}}\right)^{-1}
$$

# Modified Newton's Method 

$$
D^{k}=\left(\nabla^{2} f\left(x^{0}\right)\right)^{-1}, \quad k=0,1, \ldots
$$

provided $\nabla^{2} f\left(x^{0}\right)$ is positive definite. This method is the same as Newton's method except that to economize on overhead, the Hessian matrix is not recalculated at each iteration. A related method is obtained when the Hessian is recomputed every $p>1$ iterations.

## Discretized Newton's Method

$$
D^{k}=\left(H\left(x^{k}\right)\right)^{-1}, \quad k=0,1, \ldots
$$

where $H\left(x^{k}\right)$ is a positive definite symmetric approximation of $\nabla^{2} f\left(x^{k}\right)$, formed by using finite difference approximations of the second derivatives, based on first derivatives or values of $f$.

## Gauss-Newton Method

This method is applicable to the problem of minimizing the sum of squares of real valued functions $g_{1}, \ldots, g_{m}$, a problem often encountered in statistical data analysis and in the context of neural network training (see Section 1.5). By denoting $g=\left(g_{1}, \ldots, g_{m}\right)$, the problem is written as

$$
\begin{aligned}
& \operatorname{minimize} f(x)=\frac{1}{2}\|g(x)\|^{2}=\frac{1}{2} \sum_{i=1}^{m}\left(g_{i}(x)\right)^{2} \\
& \text { subject to } x \in \Re^{n}
\end{aligned}
$$

We choose

$$
D^{k}=\left(\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}\right)^{-1}, \quad k=0,1, \ldots
$$

assuming the matrix $\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}$ is invertible. The latter matrix is always positive semidefinite, and it is positive definite and hence invertible if and only if the matrix $\nabla g\left(x^{k}\right)$ has rank $n$ (Prop. A. 20 in Appendix A). Since

$$
\nabla f\left(x^{k}\right)=\nabla g\left(x^{k}\right) g\left(x^{k}\right)
$$

---

the Gauss-Newton method takes the form

$$
x^{k+1}=x^{k}-\alpha^{k}\left(\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}\right)^{-1} \nabla g\left(x^{k}\right) g\left(x^{k}\right)
$$

We will see in Section 1.5 that the Gauss-Newton method may be viewed as an approximation to Newton's method, particularly when the optimal value of $\|g(x)\|^{2}$ is small.

Other choices of $D^{k}$ yield the class of Quasi-Newton methods discussed in Section 1.7. There are also some interesting methods where the direction $d^{k}$ is not usually expressed as $d^{k}=-D^{k} \nabla f\left(x^{k}\right)$. Important examples are the conjugate gradient method and the coordinate descent methods discussed in Sections 1.6 and 1.8, respectively.

# Stepsize Selection 

There are a number of rules for choosing the stepsize $\alpha^{k}$ in a gradient method. We list some that are used widely in practice:

## Minimization Rule

Here $\alpha^{k}$ is such that the cost function is minimized along the direction $d^{k}$. that is, $\alpha^{k}$ satisfies

$$
f\left(x^{k}+\alpha^{k} d^{k}\right)=\min _{\alpha \geq 0} f\left(x^{k}+\alpha d^{k}\right)
$$

## Limited Minimization Rule

This is a version of the minimization rule, which is more easily implemented in many cases. A fixed scalar $s>0$ is selected and $\alpha^{k}$ is chosen to yield the greatest cost reduction over all stepsizes in the interval $[0, s]$, i.e.,

$$
f\left(x^{k}+\alpha^{k} d^{k}\right)=\min _{\alpha \in[0, s]} f\left(x^{k}+\alpha d^{k}\right)
$$

The minimization and limited minimization rules must generally be implemented with the aid of one-dimensional line search algorithms (see Appendix C). In general, the minimizing stepsize cannot be computed exactly, and in practice, the line search is stopped once a stepsize $\alpha^{k}$ satisfying some termination criterion is obtained. Some stopping criteria are discussed in Exercise 2.16 at the end of this section.

---

# Successive Stepsize Reduction - Armijo Rule 

To avoid the often considerable computation associated with the line minimization rules, it is natural to consider rules based on successive stepsize reduction. In the simplest rule of this type an initial stepsize $s$ is chosen, and if the corresponding vector $x^{k}+s d^{k}$ does not yield an improved value of $f$, that is, $f\left(x^{k}+s d^{k}\right) \geq f\left(x^{k}\right)$, the stepsize is reduced, perhaps repeatedly, by a certain factor, until the value of $f$ is improved. While this method often works in practice, it is theoretically unsound because the cost improvement obtained at each iteration may not be substantial enough to guarantee convergence to a minimum. This is illustrated in Fig. 1.2.6.

The Armijo rule is essentially the successive reduction rule just described, suitably modified to eliminate the theoretical convergence difficulty shown in Fig. 1.2.6. Here, fixed scalars $s, \beta$, and $\sigma>0$, with $\beta \in(0,1)$, and $\sigma \in(0,1)$ are chosen, and we set $\alpha^{k}=\beta^{m_{k}} s$, where $m_{k}$ is the first nonnegative integer $m$ for which

$$
f\left(x^{k}\right)-f\left(x^{k}+\beta^{m} s d^{k}\right) \geq-\sigma \beta^{m} s \nabla f\left(x^{k}\right)^{\prime} d^{k}
$$

In other words, the stepsizes $\beta^{m} s, m=0,1, \ldots$, are tried successively until the above inequality is satisfied for $m=m_{k}$. Thus, we are not satisfied with just a cost improvement: the amount of the improvement has to be sufficiently large as per the test (2.7). Figure 1.2.7 illustrates the rule.

Usually $\sigma$ is chosen close to zero, for example, $\sigma \in\left[10^{-5}, 10^{-1}\right]$. The reduction factor $\beta$ is usually chosen from $1 / 2$ to $1 / 10$ depending on the confidence we have on the quality of the initial stepsize $s$. We can always take $s=1$ and multiply the direction $d^{k}$ by a scaling factor. Many methods, such as Newton-like methods, incorporate some type of implicit scaling of the direction $d^{k}$, which makes $s=1$ a good stepsize choice (see the discussion on rate of convergence in Section 1.3). If a suitable scaling factor for $d^{k}$ is not known, one may use various ad hoc schemes to determine one. For example, a simple possibility is based on quadratic interpolation of the function

$$
g(\alpha)=f\left(x^{k}+\alpha d^{k}\right)
$$

which is the cost along the direction $d^{k}$, viewed as a function of the stepsize $\alpha$. In this scheme, we select some stepsize $\bar{\alpha}$, evaluate $g(\bar{\alpha})$, and perform the quadratic interpolation of $g$ on the basis of $g(0)=f\left(x^{k}\right), d g(0) / d \alpha=$ $\nabla f\left(x^{k}\right)^{\prime} d^{k}$, and $g(\bar{\alpha})$. If $\hat{\alpha}$ minimizes the quadratic interpolation, we replace $d^{k}$ by $\bar{d}^{k}=\bar{\alpha} d^{k}$, and we use an initial stepsize $s=1$.

## Goldstein Rule

Here, a fixed scalar $\sigma \in(0,1 / 2)$ is selected, and $\alpha^{k}$ is chosen to satisfy

$$
\sigma \leq \frac{f\left(x^{k}+\alpha^{k} d^{k}\right)-f\left(x^{k}\right)}{\alpha^{k} \nabla f\left(x^{k}\right)^{\prime} d^{k}} \leq 1-\sigma
$$

---

![[chapter_21_p26_img13.jpeg]]

Figure 1.2.6. Example of failure of the successive stepsize reduction rule for the onedimensional function

$$
f(x)= \begin{cases}\frac{3(1-x)^{2}}{4}-2(1-x), & \text { if } x>1 \\ \frac{3(1+x)^{2}}{4}-2(1+x), & \text { if } x<-1 \\ x^{2}-1, & \text { if }-1 \leq x \leq 1\end{cases}
$$

The gradient of $f$ is given by

$$
\nabla f(x)= \begin{cases}\frac{3 x}{2}+\frac{1}{2}, & \text { if } x>1 \\ \frac{3 x}{2}-\frac{1}{2}, & \text { if } x<-1 \\ 2 x, & \text { if }-1 \leq x \leq 1\end{cases}
$$

It is seen that $f$ is strictly convex, continuously differentiable, and is minimized at $x^{*}=0$. Furthermore, for any two scalars $x$, and $\tilde{x}$ we have

$$
f(x)<f(\tilde{x}) \quad \text { if and only if } \quad|x|<|\tilde{x}|
$$

We have for $x>1$

$$
x-\nabla f(x)=x-\frac{3 x}{2}-\frac{1}{2}=-\left(1+\frac{x-1}{2}\right)
$$

from which it can be verified that $|x-\nabla f(x)|<|x|$, so that $f(x-\nabla f(x))<f(x)$ and $x-\nabla f(x)<-1$. Similarly, for $x<-1$, we have $f(x-\nabla f(x))<f(x)$ and $x-\nabla f(x)>1$. Consider now the steepest descent iteration where the stepsize is successively reduced from an initial stepsize $s=1$ until descent is obtained. Let the starting point satisfy $\left|x^{0}\right|>1$. From the preceding equations, it follows that $f\left(x^{0}-\nabla f\left(x^{0}\right)\right)<f\left(x^{0}\right)$ and the stepsize $s=1$ will be accepted by the method. Thus, the next point is $x^{1}=x^{0}-\nabla f\left(x^{0}\right)$, which satisfies $\left|x^{1}\right|>1$. By repeating the preceding argument, we see that the generated sequence $\left\{x^{k}\right\}$ satisfies $\left|x^{k}\right|>1$ for all $k$, and cannot converge to the unique stationary point $x^{*}=0$. In fact, it can be shown that $\left\{x^{k}\right\}$ will have two limit points, $\bar{x}=1$ and $\bar{x}=-1$, for every $x^{0}$ with $\left|x^{0}\right|>1$.

---

![[chapter_21_p27_img14.jpeg]]

Figure 1.2.7. Line search by the Armijo rule. We start with the trial stepsize $s$ and continue with $\beta s, \beta^{2} s, \ldots$, until the first time that $\beta^{m} s$ falls within the set of stepsizes $\alpha$ satisfying the inequality

$$
f\left(x^{k}\right)-f\left(x^{k}+\alpha d^{k}\right) \geq-\sigma \alpha \nabla f\left(x^{k}\right)^{\prime} d^{k}
$$

While this set need not be an interval, it will always contain an interval of the form $[0, \delta]$ with $\delta>0$, provided $\nabla f\left(x^{k}\right)^{\prime} d^{k}<0$. For this reason the stepsize $\alpha^{k}$ chosen by the Armijo rule is well defined and will be found after a finite number of trial evaluations of $f$ at the points $\left(x^{k}+s d^{k}\right),\left(x^{k}+\beta s d^{k}\right), \ldots$
(cf. Fig. 1.2.8). It is possible to show that if $f$ is bounded below, there exists an interval of stepsizes $\alpha^{k}$ for which the relation above is satisfied. There are fairly simple algorithms for finding such a stepsize but we will not go into the details, since in practice the simpler Armijo rule seems to be universally preferred. The Goldstein rule is stated here because of its historical significance: it was the first sound proposal for a general-purpose stepsize rule that did not rely on line minimization, and it embodies the fundamental idea on which the subsequently proposed Armijo rule was based.

# Constant Stepsize 

Here a fixed stepsize $s>0$ is selected and

$$
\alpha^{k}=s, \quad k=0,1, \ldots
$$

The constant stepsize rule is very simple. However, if the stepsize is too large, divergence will occur, while if the stepsize is too small, the rate of convergence may be very slow. Thus, the constant stepsize rule is useful only for problems

---

![[chapter_21_p28_img15.jpeg]]

Figure 1.2.8. Illustration of the set of stepsizes that are acceptable in the Goldstein rule.
where an appropriate constant stepsize value is known or can be determined fairly easily.

# Diminishing Stepsize 

Here the stepsize converges to zero,

$$
\alpha^{k} \rightarrow 0
$$

This stepsize rule is different than the preceding ones in that it does not guarantee descent at each iteration, although descent becomes more likely as the stepsize diminishes. One difficulty with a diminishing stepsize is that it may become so small that substantial progress cannot be maintained, even when far from a stationary point. For this reason, we require that

$$
\sum_{k=0}^{\infty} \alpha^{k}=\infty
$$

The last condition guarantees that $\left\{x^{k}\right\}$ does not converge to a nonstationary point. Indeed, if $x^{k} \rightarrow \bar{x}$, then for any large indexes $m$ and $n(m>n)$ we have

$$
x^{m} \approx x^{n} \approx \bar{x}, \quad x^{m} \approx x^{n}-\left(\sum_{k=n}^{m-1} \alpha^{k}\right) \nabla f(\bar{x})
$$

which is a contradiction when $\bar{x}$ is nonstationary and $\sum_{k=n}^{m-1} \alpha^{k}$ can be made arbitrarily large. Generally, the diminishing stepsize rule has good theoretical

---

convergence properties under some additional assumptions that guarantee that $\alpha^{k}$ diminishes sufficiently fast to avoid oscillatory behavior of the method (see Prop. 1.2.4, and Exercises 2.13 and 2.14). The associated convergence rate tends to be slow, so this stepsize rule is used primarily in situations where slow convergence is inevitable; for example, in singular problems or when the gradient is calculated with error (see the discussion later in this section).

# Convergence Issues 

Let us now delineate the type of convergence issues that we would like to clarify. We will first discuss informally these issues and we will state and prove the more formal convergence results in Subsection 1.2.2. Given a gradient method, ideally we would like the generated sequence $\left\{x^{k}\right\}$ to converge to a global minimum. Unfortunately, however, this is too much to expect, at least when $f$ is not convex, because of the presence of local minima that are not global. Indeed a gradient method is guided downhill by the form of $f$ near the current iterate, while being oblivious of the global structure of $f$, and thus, can easily get attracted to any type of minimum, global or not. Furthermore, if a gradient method starts or lands at any stationary point, including a local maximum, it stops at that point. Thus, the most we can expect from a gradient method is that it converges to a stationary point. Such a point is a global minimum if $f$ is convex, but this need not be so for nonconvex problems. Thus, it must be recognized that gradient methods can be quite inadequate, particularly if little is known about the location and/or other properties of global minima. For such problems one must either try an often difficult and frustrating process of running a gradient method from multiple starting points, or else resort to a fundamentally different approach.

Generally, depending on the nature of the cost function $f$, the sequence $\left\{x^{k}\right\}$ generated by a gradient method need not have a limit point; in fact $\left\{x^{k}\right\}$ is typically unbounded if $f$ has no local minima. If, however, we know that the level set $\{x \mid f(x) \leq f\left(x^{0}\right)\}$ is bounded, and the stepsize is chosen to enforce descent at each iteration, then the sequence $\left\{x^{k}\right\}$ must be bounded since it belongs to this level set. It must then have at least one limit point; this is because every bounded sequence has at least one limit point (see Prop. A. 5 of Appendix A).

Even if $\left\{x^{k}\right\}$ is bounded, convergence to a single limit point may not be easy to guarantee. However, it can be shown that local minima, which are isolated stationary points (unique stationary points within some open sphere), tend to attract most types of gradient methods, that is, once a gradient method gets sufficiently close to such a local minimum, it converges to it. This is the subject of a simple and remarkably powerful result, the capture theorem, which is given in the next subsection (Prop. 1.2.5). Exercise 2.13 develops another type of convergence result for the

---

steepest descent method, which applies to the case where there are multiple nonisolated local minima. Generally, if there is a connected set of multiple global minima, it is theoretically possible for $\left\{x^{k}\right\}$ to have multiple limit points (see Exercise 2.18), but the occurance of such a phenomenon has never been documented in practice.

# Limit Points of Gradient Methods 

We now address the question of whether each limit point of a sequence $\left\{x^{k}\right\}$ generated by a gradient method is a stationary point. From the first order Taylor expansion

$$
f\left(x^{k+1}\right)=f\left(x^{k}\right)+\alpha^{k} \nabla f\left(x^{k}\right)^{\prime} d^{k}+o\left(\alpha^{k}\right)
$$

we see that if the slope of $f$ at $x^{k}$ along the direction $d^{k}$, which is $\nabla f\left(x^{k}\right)^{\prime} d^{k}$, has "substantial" magnitude, the rate of progress of the method will also tend to be substantial. If on the other hand, the directions $d^{k}$ tend to become asymptotically orthogonal to the gradient direction as $x^{k}$ approaches a nonstationary point, the slope $\nabla f\left(x^{k}\right)^{\prime} d^{k}$ will tend to zero, and there is a chance that the method will get "stuck" near that point. To ensure that this does not happen, we consider rather technical conditions on the directions $d^{k}$, which are either naturally satisfied or can be easily enforced in most algorithms of interest.

One such condition for the case where

$$
d^{k}=-D^{k} \nabla f\left(x^{k}\right)
$$

is to assume that the eigenvalues of the positive definite symmetric matrix $D^{k}$ are bounded above and bounded away from zero, that is, for some positive scalars $c_{1}$ and $c_{2}$, we have

$$
c_{1}\|z\|^{2} \leq z^{\prime} D^{k} z \leq c_{2}\|z\|^{2}, \quad \forall z \in \Re^{n}, k=0,1, \ldots
$$

It can be seen then that

$$
\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right|=\left|\nabla f\left(x^{k}\right)^{\prime} D^{k} \nabla f\left(x^{k}\right)\right| \geq c_{1}\left\|\nabla f\left(x^{k}\right)\right\|^{2}
$$

and

$$
\left\|d^{k}\right\|^{2}=\left|\nabla f\left(x^{k}\right)^{\prime}\left(D^{k}\right)^{2} \nabla f\left(x^{k}\right)\right| \leq c_{2}^{2}\left\|\nabla f\left(x^{k}\right)\right\|^{2}
$$

where we have used the fact that, from Eq. (2.8), $c_{2}$ is no less than the largest eigenvalue of $D^{k}$, and that the eigenvalues of $\left(D^{k}\right)^{2}$ are equal to the squares of the corresponding eigenvalues of $D^{k}$ (Props. A. 18 and A. 13 in Appendix A). Thus, as long as $\nabla f\left(x^{k}\right)$ does not tend to zero, $\nabla f\left(x^{k}\right)$ and $d^{k}$ cannot become asymptotically orthogonal.

We now introduce another "nonorthogonality" type of condition, which is more general than the "bounded eigenvalues" condition (2.8). Let us

---

assume that the direction $d^{k}$ is uniquely determined by the corresponding iterate $x^{k}$; that is, $d^{k}$ is obtained as a given function of $x^{k}$. We say that the direction sequence $\left\{d^{k}\right\}$ is grndient related to $\left\{x^{k}\right\}$ if the following property can be shown:

For any subsequence $\left\{x^{k}\right\}_{k \in \mathcal{K}}$ that converges to a nonstationary point, the corresponding subsequence $\left\{d^{k}\right\}_{k \in \mathcal{K}}$ is bounded and satisfies

$$
\limsup _{k \rightarrow \infty, k \in \mathcal{K}} \nabla f\left(x^{k}\right)^{\prime} d^{k}<0
$$

In particular, if $\left\{d^{k}\right\}$ is gradient related, it follows that if a subsequence $\left\{\nabla f\left(x^{k}\right)\right\}_{k \in \mathcal{K}}$ tends to a nonzero vector, the corresponding subsequence of directions $d^{k}$ is bounded and does not tend to be orthogonal to $\nabla f\left(x^{k}\right)$. Roughly, this means that $d^{k}$ does not become "too small" or "too large" relative to $\nabla f\left(x^{k}\right)$, and that the angle between $d^{k}$ and $\nabla f\left(x^{k}\right)$ does not get "too close" to 90 degrees.

We can often guarantee a priori that $\left\{d^{k}\right\}$ is gradient related. In particular, if $d^{k}=-D^{k} \nabla f\left(x^{k}\right)$ and the eigenvalues of $D^{k}$ are bounded as in the "bounded eigenvalues" condition (2.8), it can be seen that $\left\{d^{k}\right\}$ is gradient related, provided $x^{k}$ is nonstationary for all $k$ (if $x^{k}$ is stationary for some $k$, the issue of convergence in effect does not arise). Two other examples of conditions that, if satisfied for some scalars $c_{1}>0, c_{2}>0$, $p_{1} \geq 0, p_{2} \geq 0$, and all $k$, guarantee that $\left\{d^{k}\right\}$ is gradient related are
(a)

$$
c_{1}\left\|\nabla f\left(x^{k}\right)\right\|^{p_{1}} \leq-\nabla f\left(x^{k}\right)^{\prime} d^{k}, \quad\left\|d^{k}\right\| \leq c_{2}\left\|\nabla f\left(x^{k}\right)\right\|^{p_{2}}
$$

(b)

$$
d^{k}=-D^{k} \nabla f\left(x^{k}\right)
$$

with $D^{k}$ a positive definite symmetric matrix satisfying

$$
c_{1}\left\|\nabla f\left(x^{k}\right)\right\|^{p_{1}}\|z\|^{2} \leq z^{\prime} D^{k} z \leq c_{2}\left\|\nabla f\left(x^{k}\right)\right\|^{p_{2}}\|z\|^{2}, \quad \forall z \in \Re^{n}
$$

This condition generalizes the "bounded eigenvalues" condition (2.8), which is obtained for $p_{1}=p_{2}=0$.

An important convergence result is that if $\left\{d^{k}\right\}$ is gradient related and the minimization rule, or the limited minimization rule, or the Armijo rule is used, then all limit points of $\left\{x^{k}\right\}$ are stationary. This is shown in Prop. 1.2.1, given in the next subsection. Proposition 1.2.2 provides a similar result for the Goldstein rule.

When a constant stepsize is used, convergence can be proved assuming that the stepsize is sufficiently small and that $f$ satisfies some further conditions (cf. Prop. 1.2.3). Under the same conditions, convergence can also be proved for a diminishing stepsize.

---

There is a common line of proof for these convergence results. The main idea is that the cost function is improved at each iteration and that, based on our assumptions, the improvement is "substantial" near a nonstationary point, i.e., it is bounded away from zero. We then argue that the algorithm cannot approach a nonstationary point, since in this case the total cost improvement would accumulate to infinity.

# Termination of Gradient Methods 

Generally, gradient methods are not finitely convergent, so it is necessary to have criteria for terminating the iterations with some assurance that we are reasonably close to at least a local minimum. A typical approach is to stop the computation when the norm of the gradient becomes sufficiently small, that is, when a point $x^{k}$ is obtained with

$$
\left\|\nabla f\left(x^{k}\right)\right\| \leq \epsilon
$$

where $\epsilon$ is a small positive scalar. Unfortunately, it is not known a priori how small one should take $\epsilon$ in order to guarantee that the final point $x^{k}$ is a "good" approximation to a stationary point. The appropriate value of $\epsilon$ depends on how the problem is scaled. In particular, if $f$ is multiplied by some scalar, the appropriate value of $\epsilon$ is also multiplied by the same scalar. It is possible to correct this difficulty by replacing the criterion $\left\|\nabla f\left(x^{k}\right)\right\| \leq \epsilon$ with

$$
\frac{\left\|\nabla f\left(x^{k}\right)\right\|}{\left\|\nabla f\left(x^{0}\right)\right\|} \leq \epsilon
$$

Still, however, the gradient norm $\left\|\nabla f\left(x^{k}\right)\right\|$ depends on all the components of the gradient, and depending on how the optimization variables are scaled, the preceding termination criterion may not work well. In particular, some components of the gradient may be naturally much smaller than others, thus requiring a smaller value of $\epsilon$ than the other components.

Assuming that the direction $d^{k}$ captures the relative scaling of the optimization variables, it may be appropriate to terminate computation when the norm of the direction $d^{k}$ becomes sufficiently small, that is,

$$
\left\|d^{k}\right\| \leq \epsilon
$$

Still the appropriate value of $\epsilon$ may not be easy to guess, and it may be necessary to experiment prior to settling on a reasonable termination criterion for a given problem. Sometimes, other problem-dependent criteria are used, in addition or in place of the criteria $\left\|\nabla f\left(x^{k}\right)\right\| \leq \epsilon$ and $\left\|d^{k}\right\| \leq \epsilon$.

When $\nabla^{2} f(x)$ is positive definite, the condition $\left\|\nabla f\left(x^{k}\right)\right\| \leq \epsilon$ yields bounds on the distance from local minima. In particular, if $x^{*}$ is a local minimum of $f$ and there exists $m>0$ such that for all $x$ in a sphere $S$ centered at $x^{*}$ we have

$$
m\|z\|^{2} \leq z^{\prime} \nabla^{2} f(x) z, \quad \forall z \in \Re^{n}
$$

---

then every $x \in S$ satisfying $\|\nabla f(x)\| \leq \epsilon$ also satisfies

$$
\left\|x-x^{*}\right\| \leq \frac{\epsilon}{m}, \quad f(x)-f\left(x^{*}\right) \leq \frac{\epsilon^{2}}{m}
$$

(see Exercise 2.4).
In the absence of positive definiteness conditions, it may be very difficult to infer the proximity of the current iterate to the optimal solution set by just using the gradient norm. We will return to this point when we will discuss singular local minima in the next section.

# Spacer Steps 

Often, optimization problems are solved with complex descent algorithms in which the rule used to determine the next point may depend on several previous points or on the iteration index $k$. Some of the conjugate direction algorithms discussed in Section 1.6 are of this type. Other algorithms consist of a combination of different methods and switch from one method to the other in a manner that may either be prespecified or may depend on the progress of the algorithm. Such combinations are usually introduced in order to improve speed of convergence or reliability. However, their convergence analysis can become extremely complicated. It is thus often valuable to know that if in such algorithms one inserts, perhaps irregularly but infinitely often, an iteration of a convergent algorithm such as the gradient methods of this section, then the theoretical convergence properties of the overall algorithm are quite satisfactory. Such an iteration is known as a spacer step. The related convergence result is given in Prop. 1.2.6. The only requirement imposed on the iterations of the algorithm other than the spacer steps is that they do not increase the cost; these iterations, however, need not strictly decrease the cost.

## Gradient Methods with Random and Nonrandom Errors*

Frequently in optimization problems the gradient $\nabla f\left(x^{k}\right)$ is not computed exactly. Instead, one has available

$$
g^{k}=\nabla f\left(x^{k}\right)+e^{k}
$$

where $e^{k}$ is an uncontrollable error vector. There are several potential sources of error; roundoff error, and discretization error due to finite difference approximations to the gradient are two possibilities, but there are others, which will be discussed in Section 1.5. Let us for concreteness focus on the steepest descent method with errors,

$$
x^{k+1}=x^{k}-\alpha^{k} g^{k}
$$

---

and let us consider several qualitatively different cases:
(a) $e^{k}$ is small relative to the gradient, that is,

$$
\left\|e^{k}\right\|<\left\|\nabla f\left(x^{k}\right)\right\|, \quad \forall k
$$

Then, assuming $\nabla f\left(x^{k}\right) \neq 0,-g^{k}$ is a direction of cost improvement, that is, $\nabla f\left(x^{k}\right)^{\prime} g^{k}>0$. This is illustrated in Fig. 1.2.9, and is verified by the calculation

$$
\begin{aligned}
\nabla f\left(x^{k}\right)^{\prime} g^{k} & =\left\|\nabla f\left(x^{k}\right)\right\|^{2}+\nabla f\left(x^{k}\right)^{\prime} e^{k} \\
& \geq\left\|\nabla f\left(x^{k}\right)\right\|^{2}-\left\|\nabla f\left(x^{k}\right)\right\|\left\|e^{k}\right\| \\
& =\left\|\nabla f\left(x^{k}\right)\right\|\left(\left\|\nabla f\left(x^{k}\right)\right\|-\left\|e^{k}\right\|\right) \\
& >0
\end{aligned}
$$

In this case convergence results that are analogous to Props. 1.2.3 and 1.2.4 can be shown.
![[chapter_21_p34_img16.jpeg]]

Figure 1.2.9. Illustration of the descent property of the direction $g^{k}=\nabla f\left(x^{k}\right)+e^{k}$. If the error $e^{k}$ has smaller norm than the gradient $\nabla f\left(x^{k}\right)$, then $g^{k}$ lies strictly within the sphere centered at $\nabla f\left(x^{k}\right)$ with radius $\left\|\nabla f\left(x^{k}\right)\right\|$, and thus makes an angle less than 90 degrees with $\nabla f\left(x^{k}\right)$.
(b) $\left\{e^{k}\right\}$ is bounded, that is,

$$
\left\|e^{k}\right\| \leq \delta, \quad \forall k
$$

where $\delta$ is some scalar. Then by the preceding calculation (2.10), the method operates like a descent method within the region

$$
\{x \mid\|\nabla f(x)\|>\delta\}
$$

In the complementary region where $\|\nabla f(x)\| \leq \delta$, the method can behave quite unpredictably. For example, if the errors $e^{k}$ are constant, say $e^{k} \equiv e$, then since $g^{k}=\nabla f\left(x^{k}\right)+e$, the method will essentially be trying to minimize $f(x)+e^{\prime} x$ and will typically converge to a point $\bar{x}$ with $\nabla f(\bar{x})=-e$. If the errors $e^{k}$ vary substantially, the method will tend to oscillate within the region where $\|\nabla f(x)\| \leq \delta$ (see Exercise 2.17 and also Exercise 3.4 in the next section). The precise behavior will depend on the precise nature of the errors, and also on whether

---

a constant or a diminishing stepsize is used (see also the following cases).
(c) $\left\{e^{k}\right\}$ is proportional to the stepsize, that is,

$$
\left\|e^{k}\right\| \leq q \alpha^{k}, \quad \forall k
$$

where $q$ is some scalar. If the stepsize is constant, we come under case (b), while if the stepsize is diminishing, the behavior described in case (b) applies, but with $\delta \rightarrow 0$, so the method will tend to converge to a stationary point of $f$. Important situations where the condition $\left\|e^{k}\right\| \leq q \alpha^{k}$ holds will be encountered in Section 1.5 (see Prop. 1.5.2 of that section).
(d) $\left\{e^{k}\right\}$ are independent zero mean random vectors with finite variance. An important special case where such errors arise is when $f$ is of the form

$$
f(x)=E_{w}\{F(x, w)\}
$$

where $F: \Re^{m+n} \rightarrow \Re$ is some function, $w$ is a random vector in $\Re^{m}$, and $E_{w}\{\cdot\}$ denotes expected value. Under very mild assumptions it can be shown that if $F$ is continuously differentiable, the same is true of $f$ and furthermore,

$$
\nabla f(x)=E_{w}\left\{\nabla_{x} F(x, w)\right\}
$$

Often an approximation $g^{k}$ to $\nabla f\left(x^{k}\right)$ is computed by simulation or by using a limited number of samples of $\nabla F(x, w)$, with potentially substantial error resulting. In an extreme case, we have

$$
g^{k}=\nabla_{x} F\left(x^{k}, w^{k}\right)
$$

where $u^{k}$ is a single sample value corresponding to $\dot{x}^{k}$. Then the error

$$
e^{k}=\nabla_{x} F\left(x^{k}, w^{k}\right)-\nabla f\left(x^{k}\right)=\nabla_{x} F\left(x^{k}, w^{k}\right)-E_{w}\left\{\nabla_{x} F\left(x^{k}, w\right)\right\}
$$

need not diminish with $\left\|\nabla f\left(x^{k}\right)\right\|$, but has zero mean, and under appropriate conditions, its effects are "averaged out". What is happening here is that the descent condition $\nabla f\left(x^{k}\right)^{\prime} g^{k}>0$ holds on the average at nonstationary points $x^{k}$. It is still possible that for some sample values of $e^{k}$, the direction $g^{k}$ is "bad", but with a diminishing stepsize, the occasional use of a bad direction cannot deteriorate the cost enough for the method to oscillate, given that on the average the method uses "good" directions. The detailed analysis of gradient methods with random errors is beyond the scope of this text. We refer to the literature (see e.g. [BeT89], [KuC78], [PoT73a], [Pol87], [TBA86]).

---

# The Role of Convergence Analysis 

The following subsection gives a number of mathematical propositions relating to the convergence properties of gradient methods. The meaning of these propositions is usually quite intuitive but their statement often requires complicated mathematical assumptions. Furthermore, their proof often involves tedious $\epsilon-\delta$ arguments, so at first sight students may wonder whether "we really have to go through all this".

When Euclid was faced with a similar question from king Ptolemy of Alexandria, he replied that "there is no royal road to geometry". In our case, however, the answer is not so simple because we are not dealing with a pure subject such as geometry that may be developed without regard for its practical application. In the eyes of most people, the value of an analysis or algorithm in nonlinear programming is judged primarily by its practical impact in solving various types of problems. It is therefore important to give some thought to the interface between convergence analysis and its practical application. To this end it is useful to consider two extreme viewpoints; most workers in the field find themselves somewhere between the two.

In the first viewpoint, convergence analysis is considered primarily a mathematical subject. The properties of an algorithm are quantified to the extent possible through mathematical statements. General and broadly applicable assertions, and simple and elegant proofs are at a premium here. The rationale is that simple statements and proofs are more readily understood, and general statements apply not only to the problems at hand but also to other problems that are likely to appear in the future. On the negative side, one may remark that simplicity is not always compatible with relevance, and broad applicability is often achieved through assumptions that are hard to verify or appreciate.

The second viewpoint largely rejects the role of mathematical analysis. The rationale here is that the validity and the properties of an algorithm for a given class of problems must be verified through practical experimentation anyway, so if an algorithm looks promising on intuitive grounds, why bother with a convergence analysis. Furthermore, there are a number of important practical questions that are hard to address analytically, such as roundoff error, multiple local minima, and a variety of finite termination and approximation issues. The main criticism of this viewpoint is that mathematical analysis often reveals (and explains) fundamental flaws of algorithms that experimentation may miss. These flaws often point the way to better algorithms or modified algorithms that are tailored to the type of practical problem at hand. Similarly, analysis may be more effective than experimentation in delineating the types of problems for which particular algorithms are well-suited.

Our own mathematical approach is tempered by practical concerns, but we note that the balance between theory and practice in nonlinear

---

programming is particularly delicate, subjective, and problem dependent. Aside from the fact that the mathematical proofs themselves often provide valuable insight into algorithms, here are some of our reasons for insisting on a rigorous convergence analysis:

(a) We want to delineate the range of applicability of various methods. In particular, we want to know for what type of cost function (once or twice differentiable, convex or nonconvex, with singular or nonsingular minima) each algorithm is best suited. If the cost function violates the assumptions under which a given algorithm can be proved to converge, it is reasonable to suspect that the algorithm is unsuitable for this cost function.

(b) We want to provide information about the qualitative behavior of various methods. For example, we want to know whether convergence of the method depends on the availability of a good starting point, whether the iterates $x^k$ or just the function values $f(x^k)$ are guaranteed to converge, etc. This information may supplement and/or guide the computational experimentation.

(c) We want to provide guidelines for choosing a few algorithms for further experimentation out of the often bewildering array of candidate algorithms that are applicable for the solution of a given type of problem. The principal means for this is the rate of convergence analysis to be given in Section 1.3. Note here that while an algorithm may provably converge, in practice it may be entirely inappropriate for a given problem because it converges very slowly. Experience has shown that without a good understanding of the rate of convergence properties of algorithms it may be difficult to exclude bad candidates from consideration without costly experimentation.

At the same time one should be aware of some of the limitations of the mathematical results that we will provide. For example, some of the assumptions under which an algorithm will be proved convergent may be hard to verify for a given type of problem. Furthermore, our convergence rate analysis of Section 1.3 is largely asymptotic; that is, it applies near the eventual limit of the generated sequence. It is possible, that an algorithm has a good asymptotic rate of convergence but it works poorly in practice for a given type of problem because it is very slow in its initial phase.

There is still another viewpoint, which is worth addressing because it is often adopted by the casual user of nonlinear programming algorithms. This user is really interested in a particular application of nonlinear programming in his/her special field, and is counting on an existing code or package to solve the problem (several such packages are commercially or publicly available). Since the package will do most of the work, the user may hope that a superficial acquaintance with the properties of the algorithms underlying the package will suffice. This hope is sometimes realized

---

but unfortunately in the majority of the cases it is not. There are a number of reasons for this. First, there are many packages implementing a lot of different methods, and to choose the right package, one needs to have insight into the suitability of different methods for the special features of the application at hand. Second, to use a package one must often know how to suitably formulate the problem, how to set various parameters (e.g. termination criteria, stepsize parameters, etc.), and how to interpret the results of the computation (particularly when things don't work out as hoped initially, which is often the case). For this, one needs considerable insight into the inner workings of the algorithm underlying the package. Finally, for a challenging practical optimization problem (e.g. one of large dimension), it may be essential to exploit its special structure, and packages often do not have this capability. As a result the user may have to modify the package or write an altogether new code that is tailored to the application at hand. Both of these require an intimate understanding of the convergence properties and other characteristics of the relevant nonlinear programming algorithms.

# 1.2.2 Convergence Results 

We now provide an analysis of the convergence behavior of gradient methods. The following proposition is the main convergence result.

Proposition 1.2.1: (Stationarity of Limit Points for Gradient Methods) Let $\left\{x^{k}\right\}$ be a sequence generated by a gradient method $x^{k+1}=x^{k}+\alpha^{k} d^{k}$, and assume that $\left\{d^{k}\right\}$ is gradient related [cf. Eq. (2.9)] and $\alpha^{k}$ is chosen by the minimization rule, or the limited minimization rule, or the Armijo rule. Then every limit point of $\left\{x^{k}\right\}$ is a stationary point.

Proof: Consider the Armijo rule, and to arrive at a contradiction, assume that $\bar{x}$ is a limit point of $\left\{x^{k}\right\}$ with $\nabla f(\bar{x}) \neq 0$. Note that since $\left\{f\left(x^{k}\right)\right\}$ is monotonically nonincreasing, $\left\{f\left(x^{k}\right)\right\}$ either converges to a finite value or diverges to $-\infty$. Since $f$ is continuous, $f(\bar{x})$ is a limit point of $\left\{f\left(x^{k}\right)\right\}$, so it follows that the entire sequence $\left\{f\left(x^{k}\right)\right\}$ converges to $f(\bar{x})$. Hence,

$$
f\left(x^{k}\right)-f\left(x^{k+1}\right) \rightarrow 0
$$

By the definition of the Armijo rule, we have

$$
f\left(x^{k}\right)-f\left(x^{k+1}\right) \geq-\sigma \alpha^{k} \nabla f\left(x^{k}\right)^{\prime} d^{k}
$$

Hence, $\alpha^{k} \nabla f\left(x^{k}\right)^{\prime} d^{k} \rightarrow 0$. Let $\left\{x^{k}\right\}_{\mathcal{K}}$ be a subsequence converging to $\bar{x}$. Since $\left\{d^{k}\right\}$ is gradient related, we have

$$
\limsup _{\substack{k \rightarrow \infty \\ k \in \mathcal{K}}} \nabla f\left(x^{k}\right)^{\prime} d^{k}<0
$$

---

and therefore

$$
\left\{\alpha^{k}\right\}_{\mathcal{K}} \rightarrow 0
$$

Hence, by the definition of the Armijo rule, we must have for some index $\bar{k} \geq 0$

$$
f\left(x^{k}\right)-f\left(x^{k}+\left(\alpha^{k} / \beta\right) d^{k}\right)<-\sigma\left(\alpha^{k} / \beta\right) \nabla f\left(x^{k}\right)^{\prime} d^{k}, \quad \forall k \in \mathcal{K}, k \geq \bar{k}
$$

that is, the initial stepsize $s$ will be reduced at least once for all $k \in \mathcal{K}$, $k \geq \bar{k}$. Denote

$$
p^{k}=\frac{d^{k}}{\left\|d^{k}\right\|}, \quad \bar{\alpha}^{k}=\frac{\alpha^{k}\left\|d^{k}\right\|}{\beta}
$$

Since $\left\{d^{k}\right\}$ is gradient related. $\left\{\left\|d^{k}\right\|\right\}_{\mathcal{K}}$ is bounded, and it follows that

$$
\left\{\bar{\alpha}^{k}\right\}_{\mathcal{K}} \rightarrow 0
$$

Since $\left\|p^{k}\right\|=1$ for all $k \in \mathcal{K}$, there exists a subsequence $\left\{p^{k}\right\}_{\overline{\mathcal{K}}}$ of $\left\{p^{k}\right\}_{\mathcal{K}}$ such that

$$
\left\{p^{k}\right\}_{\overline{\mathcal{K}}} \rightarrow \bar{p}
$$

where $\bar{p}$ is some vector with $\|\bar{p}\|=1$ [Prop. A.5(c) in Appendix A]. From Eq. (2.13), we have

$$
\frac{f\left(x^{k}\right)-f\left(x^{k}+\bar{\alpha}^{k} p^{k}\right)}{\bar{\alpha}^{k}}<-\sigma \nabla f\left(x^{k}\right)^{\prime} p^{k}, \quad \forall k \in \overline{\mathcal{K}}, k \geq \bar{k}
$$

By using the mean value theorem, this relation is written as

$$
-\nabla f\left(x^{k}+\bar{\alpha}^{k} p^{k}\right)^{\prime} p^{k}<-\sigma \nabla f\left(x^{k}\right)^{\prime} p^{k}, \quad \forall k \in \overline{\mathcal{K}}, k \geq \bar{k}
$$

where $\bar{\alpha}^{k}$ is a scalar in the interval $\left[0, \bar{\alpha}^{k}\right]$. Taking limits in the above equation we obtain

$$
-\nabla f(\bar{x})^{\prime} \bar{p} \leq-\sigma \nabla f(\bar{x})^{\prime} \bar{p}
$$

or

$$
0 \leq(1-\sigma) \nabla f(\bar{x})^{\prime} \bar{p}
$$

Since $\sigma<1$, it follows that

$$
0 \leq \nabla f(\bar{x})^{\prime} \bar{p}
$$

On the other hand, we have

$$
\nabla f\left(x^{k}\right)^{\prime} p^{k}=\frac{\nabla f\left(x^{k}\right)^{\prime} d^{k}}{\left\|d^{k}\right\|}
$$

---

By taking the limit as $k \in \overline{\mathcal{K}}, k \rightarrow \infty$,

$$
\nabla f(\bar{x})^{\prime} \bar{p} \leq \frac{\lim \sup _{k \rightarrow \infty, k \in \overline{\mathcal{K}}} \nabla f\left(x^{k}\right)^{\prime} d^{k}}{\lim \sup _{k \rightarrow \infty, k \in \overline{\mathcal{K}}}\left\|d^{k}\right\|}<0
$$

which contradicts Eq. (2.15). This proves the result for the Armijo rule.
Consider now the minimization rule, and let $\left\{x^{k}\right\}_{\mathcal{K}}$ converge to $\bar{x}$ with $\nabla f(\bar{x}) \neq 0$. Again we have that $\left\{f\left(x^{k}\right)\right\}$ decreases monotonically to $f(\bar{x})$. Let $\tilde{x}^{k+1}$ be the point generated from $x^{k}$ via the Armijo rule, and let $\tilde{\alpha}^{k}$ be the corresponding stepsize. We have

$$
f\left(x^{k}\right)-f\left(x^{k+1}\right) \geq f\left(x^{k}\right)-f\left(\tilde{x}^{k+1}\right) \geq-\sigma \tilde{\alpha}^{k} \nabla f\left(x^{k}\right)^{\prime} d^{k}
$$

By repeating the arguments of the earlier proof following Eq. (2.12), replacing $\alpha^{k}$ by $\tilde{\alpha}^{k}$, we can obtain a contradiction. In particular, we have

$$
\left\{\tilde{\alpha}^{k}\right\}_{\mathcal{K}} \rightarrow 0
$$

and by the definition of the Armijo rule, we have for some index $\bar{k} \geq 0$

$$
f\left(x^{k}\right)-f\left(x^{k}+\left(\tilde{\alpha}^{k} / \beta\right) d^{k}\right)<-\sigma\left(\tilde{\alpha}^{k} / \beta\right) \nabla f\left(x^{k}\right)^{\prime} d^{k}, \quad \forall k \in \mathcal{K}, k \geq \bar{k}
$$

[cf. Eq. (2.13)]. Proceeding as earlier, we obtain Eqs. (2.14) and (2.15) (with $\bar{\alpha}^{k}=\tilde{\alpha}^{k}\left\|d^{k}\right\| / \beta$ ), and a contradiction of Eq. (2.15).

The line of argument just used establishes that any stepsize rule that gives a larger reduction in cost at each iteration than the Armijo rule inherits its convergence properties. This also proves the proposition for the limited minimization rule. Q.E.D.

The following proposition can be shown similar to Prop. 1.2.1. Its proof is left for the reader.

Proposition 1.2.2: The conclusions of Prop. 1.2.1 hold if $\left\{d^{k}\right\}$ is gradient related and $\alpha^{k}$ is chosen by the Goldstein rule.

The next proposition establishes, among other things, convergence for the case of a constant stepsize. The idea is that if the rate of growth of the gradient of $f$ is limited [i.e., the curvature of $f$ is limited; see the assumption (2.16) that follows], then one can construct a quadratic function $f^{k}$ that majorizes $f$; see Fig. 1.2.10. Given $x^{k}$ and $d^{k}$, an appropriate constant stepsize $\alpha^{k}$ can then be obtained within an interval around the scalar $\bar{\alpha}$ that minimizes $f^{k}$ along the direction $d^{k}$.

---

![[chapter_21_p41_img17.jpeg]]

Figure 1.2.10. The idea of the proof of Prop. 1.2.3. Given $x^{k}$ and the descent direction $d^{k}$, the cost difference $f\left(x^{k}+\alpha d^{k}\right)-f\left(x^{k}\right)$ is majorized by $\alpha \nabla f\left(x^{k}\right)^{\prime} d^{k}+$ $\frac{1}{2} \alpha^{2} L\left\|d^{k}\right\|^{2}$ [based on the Lipschitz assumption (2.16); see Eq. (2.20)]. Minimization of this function over $\alpha$ yields the stepsize

$$
\bar{\alpha}=\frac{\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right|}{L\left\|d^{k}\right\|^{2}}
$$

[cf. Eq. (2.17)]. This stepsize reduces the cost function $f$ as well [cf. Eq. (2.21)].

Proposition 1.2.3: (Convergence for a Constant Stepsize) Let $\left\{x^{k}\right\}$ be a sequence generated by a gradient method $x^{k+1}=x^{k}+\alpha^{k} d^{k}$, where $\left\{d^{k}\right\}$ is gradient related. Assume that for some constant $L>0$, we have

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|, \quad \forall x, y \in \Re^{n}
$$

and that there exists a scalar $\epsilon$ such that for all $k$ we have $d^{k} \neq 0$ and

$$
0<\epsilon \leq \alpha^{k} \leq \frac{(2-\epsilon)\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right|}{L\left\|d^{k}\right\|^{2}}
$$

Then every limit point of $\left\{x^{k}\right\}$ is a stationary point of $f$.

Note: If $\left\{d^{k}\right\}$ is such that there exist positive scalars $c_{1}, c_{2}$ such that for all $k$ we have

$$
c_{1}\left\|\nabla f\left(x^{k}\right)\right\|^{2} \leq-\nabla f\left(x^{k}\right)^{\prime} d^{k}, \quad\left\|d^{k}\right\|^{2} \leq c_{2}\left\|\nabla f\left(x^{k}\right)\right\|^{2}
$$

then Eq. (2.17) is satisfied if for all $k$ we have

$$
0<\epsilon \leq \alpha^{k} \leq \frac{c_{1}(2-\epsilon)}{L c_{2}}
$$

---

In particular, for steepest descent $\left[d^{k}=-\nabla f\left(x^{k}\right)\right]$, we can take $c_{1}=c_{2}=1$, and the condition on the stepsize becomes

$$
0<\epsilon \leq \alpha^{k} \leq \frac{2-\epsilon}{L}
$$

Furthermore, if for all $k$, we have $d^{k}=-D^{k} \nabla f\left(x^{k}\right)$ with $D^{k}$ positive definite symmetric and having eigenvalues in an interval $[\gamma, \Gamma]$, the condition (2.18) can be seen to hold with

$$
c_{1}=\gamma, \quad c_{2}=\Gamma^{2}
$$

Proof: By using the descent lemma (Prop. A. 24 of Appendix A), we obtain

$$
\begin{aligned}
f\left(x^{k}+\alpha^{k} d^{k}\right)-f\left(x^{k}\right) & \leq \alpha^{k} \nabla f\left(x^{k}\right)^{\prime} d^{k}+\frac{1}{2}\left(\alpha^{k}\right)^{2} L\left\|d^{k}\right\|^{2} \\
& =\alpha^{k}\left(\frac{1}{2} \alpha^{k} L\left\|d^{k}\right\|^{2}-\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right|\right)
\end{aligned}
$$

The right-hand side of Eq. (2.17) yields

$$
\frac{1}{2} \alpha^{k} L\left\|d^{k}\right\|^{2}-\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right| \leq-\frac{1}{2} \epsilon\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right|
$$

Using this relation together with the condition $\alpha^{k} \geq \epsilon$ in the inequality (2.20), we obtain the following bound on the cost improvement obtained at iteration $k$ :

$$
f\left(x^{k}\right)-f\left(x^{k}+\alpha^{k} d^{k}\right) \geq \frac{1}{2} \epsilon^{2}\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right|
$$

Now if a subsequence $\left\{x^{k}\right\}_{\mathcal{K}}$ converges to a nonstationary point $\bar{x}$, we must have, as in the proof of Prop. 1.2.1, $f\left(x^{k}\right)-f\left(x^{k+1}\right) \rightarrow 0$, and Eq. (2.21) implies that $\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right| \rightarrow 0$. This contradicts the hypothesis that $\left\{d^{k}\right\}$ is gradient related. Hence, every limit point of $\left\{x^{k}\right\}$ is stationary. Q.E.D.

A condition of the form

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|, \quad \forall x, y \in \Re^{n}
$$

[cf. Eq. (2.16)] is called a Lipschitz continuity condition on $\nabla f$, and requires roughly that the "curvature" of $f$ is no more than $L$ in all directions. In particular, it is possible to show that this condition is satisfied for some $L>0$, if $f$ is twice differentiable and the Hessian $\nabla^{2} f$ is bounded over $\Re^{n}$. Unfortunately, however, it is generally difficult to obtain an estimate of $L$, so in most cases the interval of stepsizes that guarantee convergence [cf. Eq. (2.17) or (2.19)] is unknown. Thus, experimentation may be necessary to obtain an appropriate range of stepsizes.

Exercise 2.3 provides an example showing that the Lipschitz continuity condition (2.16) is essential for the validity of Prop. 1.2.3. Exercise 2.5

---

shows that this condition can be weakened somewhat; it is sufficient that it holds for all $x, y$ in the set $\left\{z \mid f(z) \leq f\left(x^{0}\right)\right\}$. The Lipschitz continuity condition also essentially guarantees convergence for a diminishing stepsize, as shown by the following proposition.

Proposition 1.2.4: (Convergence for a Diminishing Stepsize) Let $\left\{x^{k}\right\}$ be a sequence generated by a gradient method $x^{k+1}=x^{k}+$ $\alpha^{k} d^{k}$. Assume that for some constant $L>0$, we have

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|, \quad \forall x, y \in \Re^{n}
$$

and that there exist positive scalars $c_{1}, c_{2}$ such that for all $k$ we have

$$
c_{1}\left\|\nabla f\left(x^{k}\right)\right\|^{2} \leq-\nabla f\left(x^{k}\right)^{\prime} d^{k}, \quad\left\|d^{k}\right\|^{2} \leq c_{2}\left\|\nabla f\left(x^{k}\right)\right\|^{2}
$$

Suppose also that

$$
\alpha^{k} \rightarrow 0, \quad \sum_{k=0}^{\infty} \alpha^{k}=\infty
$$

Then either $f\left(x^{k}\right) \rightarrow-\infty$ or else $\left\{f\left(x^{k}\right)\right\}$ converges to a finite value and $\liminf _{k \rightarrow \infty}\left\|\nabla f\left(x^{k}\right)\right\|=0$. Furthermore, if the level set $\{x \mid$ $f(x) \leq \gamma\}$ is bounded for each scalar $\gamma$, then $\left\{x^{k}\right\}$ is bounded and every limit point of $\left\{x^{k}\right\}$ is a stationary point of $f$.

Proof: Combining Eqs. (2.20) and (2.22), we have

$$
\begin{aligned}
f\left(x^{k+1}\right) & \leq f\left(x^{k}\right)+\alpha^{k}\left(\frac{1}{2} \alpha^{k} L\left\|d^{k}\right\|^{2}-\left|\nabla f\left(x_{,}^{k}\right)^{\prime} d^{k}\right|\right) \\
& \leq f\left(x^{k}\right)-\alpha^{k}\left(c_{1}-\frac{1}{2} \alpha^{k} c_{2} L\right)\left\|\nabla f\left(x^{k}\right)\right\|^{2}
\end{aligned}
$$

Since the linear term in $\alpha^{k}$ dominates the quadratic term in $\alpha^{k}$ for sufficiently small $\alpha^{k}$, and $\alpha^{k} \rightarrow 0$, we have for some positive constant $c$ and all $k$ greater than some index $\bar{k}$,

$$
f\left(x^{k+1}\right) \leq f\left(x^{k}\right)-\alpha^{k} c\left\|\nabla f\left(x^{k}\right)\right\|^{2}
$$

From this relation, we see that for $k \geq \bar{k},\left\{f\left(x^{k}\right)\right\}$ is monotonically decreasing, so either $f\left(x^{k}\right) \rightarrow-\infty$ or $\left\{f\left(x^{k}\right)\right\}$ converges to a finite value. In the latter case, by adding Eq. (2.23) over all $k \geq \bar{k}$, we obtain

$$
c \sum_{k=\bar{k}}^{\infty} \alpha^{k}\left\|\nabla f\left(x^{k}\right)\right\|^{2} \leq f\left(x^{\bar{k}}\right)-\lim _{k \rightarrow \infty} f\left(x^{k}\right)<\infty
$$

---

We see that there cannot exist an $\epsilon>0$ such that $\left\|\nabla f\left(x^{k}\right)\right\|^{2}>\epsilon$ for all $k$, since this would contradict the assumption $\sum_{k=0}^{\infty} \alpha^{k}=\infty$. Therefore, we must have $\liminf _{k \rightarrow \infty}\left\|\nabla f\left(x^{k}\right)\right\|=0$.

Assume that the level sets of $f$ are all bounded. Since $\left\{f\left(x^{k}\right)\right\}$ is monotonically decreasing for $k \geq \bar{k}$, the sequence $\left\{x^{k}\right\}$ must belong to one of the level sets of $f$, so $\left\{x^{k}\right\}$ must be bounded. Therefore, $\left\{f\left(x^{k}\right)\right\}$ must converge to a finite value and $\liminf _{k \rightarrow \infty}\left\|\nabla f\left(x^{k}\right)\right\|=0$. We will show that $\left\|\nabla f\left(x^{k}\right)\right\| \rightarrow 0$

Indeed, assume the contrary; that is, there exists an $\epsilon>0$ such that $\left\|\nabla f\left(x^{k}\right)\right\|>\epsilon$ for all $k$ in an infinite subset of integers $\mathcal{K}$. For each $k \in \mathcal{K}$, let $i(k)$ be the first index $i$ such that $i>k$ and $\left\|\nabla f\left(x^{i}\right)\right\|<\epsilon / 2$, so that

$$
\begin{aligned}
\frac{\epsilon}{2} & \leq\left\|\nabla f\left(x^{k}\right)\right\|-\left\|\nabla f\left(x^{i(k)}\right)\right\| \leq\left\|\nabla f\left(x^{k}\right)-\nabla f\left(x^{i(k)}\right)\right\| \\
& \leq L\left\|x^{k}-x^{i(k)}\right\| \leq L \sum_{i=k}^{i(k)-1} \alpha^{i}\left\|d^{i}\right\|
\end{aligned}
$$

Since $\left\{x^{k}\right\}$ is bounded, $\left\{\nabla f\left(x^{k}\right)\right\}$ is also bounded, and the hypothesis $\left\|d^{k}\right\|^{2} \leq c_{2}\left\|\nabla f\left(x^{k}\right)\right\|^{2}$ implies that $\left\{d^{k}\right\}$ is bounded as well. Therefore, Eq. (2.24) implies that for some constant $b>0$,

$$
\frac{\epsilon}{2} \leq b \sum_{i=k}^{i(k)-1} \alpha^{i}, \quad \forall k \in \mathcal{K}
$$

Using Eq. (2.23) for sufficiently large $k \in \mathcal{K}$, and the relation $\left\|\nabla f\left(x^{i}\right)\right\| \geq$ $\epsilon / 2$ for $i=k, k+1, \ldots, i(k)-1$, we have

$$
f\left(x^{i(k)}\right) \leq f\left(x^{k}\right)-c\left(\frac{\epsilon}{2}\right)^{2} \sum_{i=k}^{i(k)-1} \alpha^{i}, \quad \forall k \in \mathcal{K}
$$

Since $\left\{f\left(x^{k}\right)\right\}$ converges to a finite value, it follows that

$$
\lim _{k \rightarrow \infty, k \in \mathcal{K}} \sum_{i=k}^{i(k)-1} \alpha^{i}=0
$$

contradicting the earlier conclusion $\frac{\epsilon}{2} \leq b \sum_{i=k}^{i(k)-1} \alpha^{i}$ for all $k \in \mathcal{K}$. Therefore, $\left\|\nabla f\left(x^{k}\right)\right\| \rightarrow 0$, which implies that every limit point of $\left\{x^{k}\right\}$ is a stationary point of $f$. Q.E.D.

Under the assumptions of the preceding proposition, descent is not guaranteed in the initial iterations. However, if the stepsizes are all sufficiently small [e.g., they satisfy the right-hand side inequality of Eq. (2.19)], descent is guaranteed at all iterations. In this case, it is sufficient that

---

the Lipschitz condition $\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|$ holds for all $x, y$ in the set $\{z \mid f(z) \leq f\left(x^{0}\right)\}$ (see Exercise 2.5): otherwise the Lipschitz condition must hold over a set larger than $\{z \mid f(z) \leq f\left(x^{0}\right)\}$ to guarantee convergence (see Exercise 2.15).

The following proposition explains to some extent why sequences generated by gradient methods tend in practice to have unique limit points. It states that local minima which happen to be isolated stationary points tend to attract gradient methods: once the method gets close enough to such a minimum it remains close and converges to it.

Proposition 1.2.5: (Capture Theorem) Let $f$ be continuously differentiable and let $\left\{x^{k}\right\}$ be a sequence satisfying $f\left(x^{k+1}\right) \leq f\left(x^{k}\right)$ for all $k$ and generated by a gradient method $x^{k+1}=x^{k}+\alpha^{k} d^{k}$, which is convergent in the sense that every limit point of sequences that it generates is a stationary point of $f$. Assume that there exist scalars $s>0$ and $c>0$ such that for all $k$ there holds

$$
\alpha^{k} \leq s, \quad\left\|d^{k}\right\| \leq c\left\|\nabla f\left(x^{k}\right)\right\|
$$

Let $x^{*}$ be a local minimum of $f$, which is the only stationary point of $f$ within some open set. Then there exists an open set $S$ containing $x^{*}$ such that if $x^{\bar{k}} \in S$ for some $\bar{k} \geq 0$, then $x^{k} \in S$ for all $k \geq \bar{k}$ and $\left\{x^{k}\right\} \rightarrow x^{*}$. Furthermore, given any scalar $\bar{\epsilon}>0$, the set $S$ can be chosen so that $\left\|x-x^{*}\right\|<\bar{\epsilon}$ for all $x \in S$.

Note: The conditions $f\left(x^{k+1}\right) \leq f\left(x^{k}\right)$ and $\alpha^{k} \leq s$ are satisfied for the Armijo rule and the limited minimization rule. They are also satisfied for a constant and a diminishing stepsize under conditions that guarantee descent at each iteration (see the proofs of Props. 1.2.3 and 1.2.4). The condition $\left\|d^{k}\right\| \leq c\|\nabla f\left(x^{k}\right)\|$ is satisfied if $d^{k}=-D^{k} \nabla f\left(x^{k}\right)$ with the eigenvalues of $D^{k}$ bounded from above.
Proof: Suppose that $\rho>0$ is such that

$$
f\left(x^{*}\right)<f(x), \quad \forall x \text { with }\left\|x-x^{*}\right\| \leq \rho
$$

Define for $t \in[0, \rho]$

$$
\phi(t)=\min _{\{x \mid t \leq\left\|x-x^{*}\right\| \leq \rho\}} f(x)-f\left(x^{*}\right)
$$

and note that $\phi$ is a monotonically nondecreasing function of $t$, and that $\phi(t)>0$ for all $t \in(0, \rho]$. Given any $\epsilon \in(0, \rho]$, let $r \in(0, \epsilon]$ be such that

$$
\left\|x-x^{*}\right\|<r \quad \Rightarrow \quad\left\|x-x^{*}\right\|+s c\|\nabla f(x)\|<\epsilon
$$

---

Consider the open set

$$
S=\left\{x \mid\left\|x-x^{*}\right\|<\epsilon, f(x)<f\left(x^{*}\right)+\phi(r)\right\}
$$

We claim that if $x^{k} \in S$ for some $k$, then $x^{k+1} \in S$.
Indeed if $x^{k} \in S$, from the definition of $\phi$ and $S$ we have

$$
\phi\left(\left\|x^{k}-x^{*}\right\|\right) \leq f\left(x^{k}\right)-f\left(x^{*}\right)<\phi(r)
$$

Since $\phi$ is monotonically nondecreasing, the above relation implies that $\left\|x^{k}-x^{*}\right\|<r$, so that by Eq. (2.25),

$$
\left\|x^{k}-x^{*}\right\|+s c\left\|\nabla f\left(x^{k}\right)\right\|<\epsilon
$$

We also have by using the hypotheses $\alpha^{k} \leq s$ and $\left\|d^{k}\right\| \leq c\left\|\nabla f\left(x^{k}\right)\right\|$

$$
\left\|x^{k+1}-x^{*}\right\| \leq\left\|x^{k}-x^{*}\right\|+\left\|\alpha^{k} d^{k}\right\| \leq\left\|x^{k}-x^{*}\right\|+s c\left\|\nabla f\left(x^{k}\right)\right\|
$$

so from the last two relations it follows that $\left\|x^{k+1}-x^{*}\right\|<\epsilon$. Since $f\left(x^{k+1}\right)<f\left(x^{k}\right)$, we also obtain $f\left(x^{k+1}\right)-f\left(x^{*}\right)<\phi(r)$, so we conclude that $x^{k+1} \in S$.

By using induction it follows that if $x^{\bar{k}} \in S$ for some $\bar{k}$, we have $x^{k} \in S$ for all $k \geq \bar{k}$. Let $\bar{S}$ be the closure of $S$. Since $\bar{S}$ is compact, the sequence $\left\{x^{k}\right\}$ will have at least one limit point, which by assumption must be a stationary point of $f$. Now the only stationary point of $f$ within $\bar{S}$ is the point $x^{*}$ (since we have $\left\|x-x^{*}\right\| \leq \rho$ for all $x \in \bar{S}$ ). Hence $x^{k} \rightarrow x^{*}$. Finally given any $\bar{\epsilon}>0$, we can choose $\epsilon \leq \bar{\epsilon}$ in which case we have $\left\|x-x^{*}\right\|<\bar{\epsilon}$ for all $x \in S$. Q.E.D.

The proof of the following proposition is similar to the one of Prop. 1.2.1, and is left for the reader.

Proposition 1.2.6: (Convergence for Spacer Steps) Consider a sequence $\left\{x^{k}\right\}$ such that

$$
f\left(x^{k+1}\right) \leq f\left(x^{k}\right), \quad k=0,1, \ldots
$$

Assume that there exists an infinite set $\mathcal{K}$ of integers for which

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}, \quad \forall k \in \mathcal{K}
$$

where $\left\{d^{k}\right\}_{\mathcal{K}}$ is gradient related and $\alpha^{k}$ is chosen by the minimization rule, or the limited minimization rule, or the Armijo rule. Then every limit point of the subsequence $\left\{x^{k}\right\}_{\mathcal{K}}$ is a stationary point.

---

# EXERCISES 

## 2.1

Consider the problem of minimizing the function of two variables $f(x, y)=$ $3 x^{2}+y^{4}$.
(a) Apply one iteration of the steepest descent method with $(1,-2)$ as the starting point and with the stepsize chosen by the Armijo rule with $s=1, \sigma=0.1$, and $d=0.5$.
(b) Repeat (a) using $s=1, \sigma=0.1, \beta=0.1$ instead. How does the cost of the new iterate compare to that obtained in (a)? Comment on the tradeoffs involved in the choice of $\beta$.
(c) Apply one iteration of Newton's method with the same starting point and stepsize rule as in (a). How does the cost of the new iterate compare to that obtained in (a)? How about the amount of work involved in finding the new iterate?

## 2.2

Describe the behavior of the steepest descent method with constant stepsize $s$ for the function $f(x)=\|x\|^{2+\beta}$, where $\beta \geq 0$. For which values of $s$ and $x^{0}$ does the method converge to $x^{*}=0$. Relate your answer to the assumptions of Prop. 1.2.3.

## 2.3

Consider the function $f: \Re^{n} \mapsto \Re$ given by

$$
f(x)=\|x\|^{3 / 2}
$$

and the method of steepest descent with a constant stepsize. Show that for this function, the Lipschitz condition $\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|$ for all $x$ and $y$ is not satisfied for any $L$. Furthermore, for any value of constant stepsize, the method either converges in a finste number of iterations to the minimizing point $x^{*}=0$ or else it does not converge to $x^{*}$.

## 2.4

Let $f$ be twice continuously differentiable. Suppose that $x^{*}$ is a local minimum such that for all $x$ in an open sphere $S$ centered at $x^{*}$, we have, for some $m>0$,

$$
m\|d\|^{2} \leq d^{\prime} \nabla^{2} f(x) d, \quad \forall d \in \Re^{n}
$$

---

Show that for every $x \in S$, we have

$$
\left\|x-x^{*}\right\| \leq \frac{\|\nabla f(x)\|}{m}, \quad f(x)-f\left(x^{*}\right) \leq \frac{\|\nabla f(x)\|^{2}}{m}
$$

Hint: Use the relation

$$
\nabla f(y)=\nabla f(x)+\int_{0}^{1} \nabla^{2} f(x+t(y-x))(y-x) d t
$$

2.5

Suppose that the Lipschitz condition

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|, \quad \forall x, y \in \Re^{n}
$$

[cf. Eq. (2.16)] is replaced by the following two conditions:
(i) For every bounded set $A \subset \Re^{n}$, there exists some constant $L$ such that

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|, \quad \forall x, y \in A
$$

(ii) The set $\{x \mid f(x) \leq c\}$ is bounded for every $c \in \Re$.

Show that:
(a) Condition (i) is always satisfied if $f$ is twice continuously differentiable.
(b) The convergence result of Prop. 1.2.3 remains valid provided that the constant stepsize $s$ is allowed to depend on the choice of the initial vector $x^{0}$. Hint: Choose a stepsize that guarantees that $x^{k}$ stays within the level set $\left\{x \mid f(x) \leq f\left(x^{0}\right)\right\}$ for all $k$.

# 2.6 

Suppose that $f$ is quadratic and of the form $f(x)=\frac{1}{2} x^{\prime} Q x-b^{\prime} x$, where $Q$ is positive definite and symmetric.
(a) Show that the Lipschitz condition $\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|$ is satisfied with $L$ equal to the maximal eigenvalue of $Q$.
(b) Consider the gradient method $x^{k+1}=x^{k}-s D \nabla f\left(x^{k}\right)$, where $D$ is positive definite and symmetric. Show that the method converges to $x^{*}=Q^{-1} b$ if $s \in(0,2 / \bar{L})$, where $\bar{L}$ is the maximum eigenvalue of $D^{1 / 2} Q D^{1 / 2}$.

---

# 2.7 

Apply the steepest descent method with constant stepsize $\alpha$ to the function $f$ of Exercise 1.11 of Section 1.1. Show that the gradient $\nabla f$ satisfies the Lipschitz condition

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|, \quad \forall x, y \in \Re
$$

for some constant $L$. Write a computer program to verify that the method is a descent method for $\alpha \in(0,2 / L)$. Do you expect to get in the limit the global minimum $x^{*}=0$ ?

## 2.8

Consider the gradient method $x^{k+1}=x^{k}+\alpha^{k} d^{k}$, where $\alpha^{k}$ is chosen by the Armijo rule or the line minimization rule and

$$
d^{k}=-\left[\begin{array}{c}
0 \\
\vdots \\
0 \\
\frac{\partial f\left(x^{k}\right)}{\partial x_{i}} \\
\vdots \\
0
\end{array}\right]
$$

where $i$ is the index for which $\left|\partial f\left(x^{k}\right) / \partial x_{j}\right|$ is maximized over $j=1, \ldots, n$. Show that every limit point of $\left\{x^{k}\right\}$ is stationary.

## 2.9

Consider the gradient method $x^{k+1}=x^{k}+\alpha^{k} d^{k}$ for the case where $f$ is positive definite quadratic, and let $\bar{\alpha}^{k}$ be the stepsize corresponding to the line minimization rule. Show that a stepsize $\alpha^{k}$ satisfies the inequalities of the Goldstein rule if and only if

$$
2 \sigma \alpha^{k} \leq \alpha^{k} \leq 2(1-\sigma) \alpha^{k}
$$

### 2.10 (Alternative Assumptions for Convergence)

Consider the gradient method $x^{k+1}=x^{k}+\alpha^{k} d^{k}$. Instead of $\left\{d^{k}\right\}$ being gradient related, assume one of the following two conditions:
(1) It can be shown that for any subsequence $\left\{x^{k}\right\}_{k \in \mathcal{K}}$ that converges to a nonstationary point, the corresponding subsequence $\left\{d^{k}\right\}_{k \in \mathcal{K}}$ is bounded and satisfies

$$
\liminf _{k \rightarrow \infty, k \in \mathcal{K}} \nabla f\left(x^{k}\right)^{\prime} d^{k}<0
$$

---

(2) $\alpha^{k}$ is chosen by the minimization rule, and for some $c>0$ and all $k$, we have

$$
\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right| \geq c\left\|\nabla f\left(x^{k}\right)\right\|\left\|d^{k}\right\| .
$$

Show that the result of Prop. 1.2.1 holds.

# 2.11 

An electrical engineer wants to maximize the current $I$ between two points A and B of a complex network by adjusting the values $x_{1}$ and $x_{2}$ of two variable resistors, where $0 \leq x_{1} \leq R_{1}, 0 \leq x_{2} \leq R_{2}$, and $R_{1}, R_{2}$ are given. The engineer does not have an adequate mathematical model of the network and decides to adopt the following procedure. She keeps the value $x_{2}$ of the second resistor fixed and adjusts the value of the first resistor until the current $I$ is maximized. She then keeps the value $x_{1}$ of the first resistor fixed and adjusts the value of the second resistor until the current $I$ is maximized. She then repeats the procedure until no further progress can be made. She knows a priori that during this procedure, the values $x_{1}$ and $x_{2}$ can never reach their extreme values $0, R_{1}$, and $R_{2}$. Explain whether there is a sound theoretical basis for the engineer's procedure. Hint: Consider how the steepest descent method works for two-dimensional problems.

### 2.12 (Behavior of Steepest Descent Near a Saddle Point)

Let $f=(1 / 2) x^{\prime} Q x$, where $Q$ is invertible and has at least one negative eigenvalue. Consider the steepest descent method with constant stepsize and show that unless the starting point $x^{0}$ belongs to the subspace spanned by the eigenvectors of $Q$ corresponding to the nonnegative eigenvalues, the generated sequence $\left\{x^{k}\right\}$ diverges.

### 2.13 (Convergence to a Single Limit)

Consider the steepest descent method $x^{k+1}=x^{k}-\alpha^{k} \nabla f\left(x^{k}\right)$ and assume that for all $x, y$, we have

$$
\frac{\|\nabla f(x)-\nabla f(y)\|^{2}}{L} \leq(\nabla f(x)-\nabla f(y))^{\prime}(x-y)
$$

(It can be shown that this condition holds if $f$ is convex, twice continuously differentiable and its Hessian matrix has eigenvalues that are less than or equal to $L$.) Assume also that $f$ has at least one stationary point. Show that $\left\{x^{k}\right\}$ converges to a stationary point of $f$ under one of the following two conditions:
(i) For some $\epsilon>0$, we have

$$
\epsilon \leq \alpha^{k} \leq \frac{2-\epsilon}{L}, \quad \forall k
$$

---

(ii) $\alpha^{k} \rightarrow 0$ and $\sum_{k=0}^{\infty} \alpha^{k}=\infty$.

Hint: Show that for any stationary point $\bar{x}$ we have

$$
\left\|x^{k+1}-\bar{x}\right\|^{2} \leq\left\|x^{k}-\bar{x}\right\|^{2}-\alpha^{k}\left(\frac{2}{L}-\alpha^{k}\right)\left\|\nabla f\left(x^{k}\right)\right\|^{2}
$$

# 2.14 (Steepest Descent with Diminishing Stepsize [CoL94]) 

Consider the steepest descent method

$$
x^{k+1}=x^{k}-\alpha^{k} \nabla f\left(x^{k}\right)
$$

assuming that the function $f$ is convex.
(a) Use the convexity of $f$ to show that for any $y \in \Re^{n}$, we have

$$
\left\|x^{k+1}-y\right\|^{2} \leq\left\|x^{k}-y\right\|^{2}-2 \alpha^{k}\left(f\left(x^{k}\right)-f(y)\right)+\left(\alpha^{k}\left\|\nabla f\left(x^{k}\right)\right\|\right)^{2}
$$

(b) Assume that

$$
\sum_{k=0}^{\infty} \alpha^{k}=\infty, \quad \alpha^{k}\left\|\nabla f\left(x^{k}\right)\right\|^{2} \rightarrow 0
$$

Show that $\liminf _{k \rightarrow \infty} f\left(x^{k}\right)=\inf _{x \in \Re^{n}} f(x)$. Hint: Argue by contradiction. Assume that for some $\delta>0$, there exists $y$ with $f(y)<f\left(x^{k}\right)-\delta$ for all $k$ sufficiently large. Use part (a).
(c) Assume that

$$
\alpha^{k}=\frac{s^{k}}{\left\|\nabla f\left(x^{k}\right)\right\|}
$$

where

$$
\sum_{k=0}^{\infty} s^{k}=\infty, \quad \sum_{k=0}^{\infty}\left(s^{k}\right)^{2}<\infty
$$

Show that $\liminf _{k \rightarrow \infty} f\left(x^{k}\right)=\inf _{x \in \Re^{n}} f(x)$, and that if $f$ has at least one global minimum, then $\left\{x^{k}\right\}$ converges to some global minimum. Hint: In part (a), set $y$ to some $x^{*}$ such that $f\left(x^{*}\right)<f\left(x^{k}\right)$ for all $k$ (if no such $x^{*}$ exists, we are done). Show that the relation

$$
\left\|x^{k+1}-x^{*}\right\|^{2} \leq\left\|x^{k}-x^{*}\right\|^{2}+\left(s^{k}\right)^{2}
$$

implies that $\left\{x^{k}\right\}$ is bounded and hence also that $\left\{\nabla f\left(x^{k}\right)\right\}$ is bounded. Use part (b).

---

# 2.15 (Divergence with Diminishing Stepsize) 

Consider the one-dimensional function

$$
f(x)=\frac{2}{3}|x|^{3}+\frac{1}{2} x^{2}
$$

and the method of steepest descent with stepsize $\alpha^{k}=\gamma /(k+1)$, where $\gamma$ is a positive scalar.
(a) Show that for $\gamma=1$ and $\left|x^{0}\right| \geq 1$ the method diverges. In particular, show that $\left|x^{k}\right| \geq k+1$ for all $k$.
(b) Characterize as best as you can the set of pairs $\left(\gamma, x^{0}\right)$ for which the method converges to $x^{*}=0$.
(c) How do you reconcile the results of (a) and (b) with Prop. 1.2.4.

### 2.16

There are several criteria for implementing approximately the minimization rule in a gradient method. An example of such a criterion is that $\alpha^{k}$ satisfies simultaneously

$$
\begin{aligned}
f\left(x^{k}\right)-f\left(x^{k}+\alpha^{k} d^{k}\right) & \geq-\sigma \alpha^{k} \nabla \bar{f}\left(x^{k}\right)^{\prime} d^{k} \\
\nabla f\left(x^{k}+\alpha^{k} d^{k}\right)^{\prime} d^{k} & \geq \beta \nabla f\left(x^{k}\right)^{\prime} d^{k}
\end{aligned}
$$

where $\alpha$ and $\beta$ are some scalars with $\sigma \in(0,1 / 2)$ and $\beta \in(\sigma, 1)$. If $\alpha^{k}$ is indeed a minimizing stepsize, then $\nabla f\left(x^{k}+\alpha^{k} d^{k}\right)^{\prime} d^{k}=d g\left(\alpha^{k}\right) / d \alpha=0$, where $g$ is the function $g(\alpha)=f\left(x^{k}+\alpha d^{k}\right)$, so Eq. (2.27) is in effect a test on the accuracy of the minimization (see Fig. 1.2.11).
(a) Show that if conditions (2.26) and (2.27) are satisfied by a gradient method at each iteration and the direction sequence is gradient related, then all limit points of the generated sequence $\left\{x^{k}\right\}$ are stationary points of $f$.
(b) Assume that there is a scalar $M$ such that $f(x) \geq M$. Show that there exists an interval $\left[c_{1}, c_{2}\right]$ with $0<c_{1}<c_{2}$, such that every $\alpha \in\left[c_{1}, c_{2}\right]$ satisfies Eqs. (2.26) and (2.27).

### 2.17 (Steepest Descent with Errors)

Consider the steepest descent method $x^{k+1}=x^{k}-\alpha^{k}\left(\nabla f\left(x^{k}\right)+e^{k}\right)$, where $e^{k}$ is an error satisfying $\left\|e^{k}\right\| \leq \delta$ for all $k$. Show that for any $\delta^{\prime}>\delta$, there exists a range of positive stepsizes $[\underline{\alpha}, \bar{\alpha}]$ such that if $\alpha^{k} \in[\underline{\alpha}, \bar{\alpha}]$ for all $k$, then either $f\left(x^{k}\right) \rightarrow-\infty$ or there is at least one limit point of $\left\{x^{k}\right\}$ in the set $\left\{x \mid\left\|\nabla f\left(x^{k}\right)\right\|<\delta^{\prime}\right\}$. Hint: Use Prop. 1.2.3.

---

![[chapter_21_p53_img18.jpeg]]

Figure 1.2.11. Illustration of the stepsize selection criterion based on Eqs. (2.26) and (2.27).

# 2.18 (A Continuum of Limit Points for Steepest Descent [Zou76]) 

Consider the two-dimensional function

$$
f(x)= \begin{cases}(r-1)^{2}-\frac{1}{2}(r-1)^{2} \cos \left(\frac{1}{r-1}-\phi\right) & \text { if } r \neq 1 \\ 0 & \text { if } r=1\end{cases}
$$

where

$$
r=\sqrt{x_{1}^{2}+x_{2}^{2}}, \quad \phi=\arctan \left(x_{1} / x_{2}\right)
$$

This function is minimized at each point of the circle where $r=1$. Consider a nonoptimal starting point and the method of steepest descent where $x^{k+1}$ is set equal to the first local minimum along the line $\left\{x^{k}-\alpha \nabla f\left(x^{k}\right) \mid \alpha \geq 0\right\}$. Show that this method follows a spiral path that comes arbitrarily close to every point of the circle of optimal points.

### 2.19 (Simplified Steepest Descent)

(a) Consider the unconstrained minimization of a function $f$ of the form

$$
f(x)=F(x, g(x))
$$

where $g: \Re^{m} \rightarrow \Re^{n}$ is continuously differentiable and $F(x, y)$ is a continuously differentiable function of the two arguments $x \in \Re^{n}$ and $y \in \Re^{m}$. It is sometimes convenient to approximate the gradient of $F(x, g(x))$ by neglecting the dependence on $g$. This leads to the method

$$
x^{k+1}=x^{k}-\alpha^{k} \nabla_{x} F\left(x^{k}, g\left(x^{k}\right)\right)
$$

---

where $\alpha^{k}$ is chosen by the minimization rule or the Armijo rule on the function $f$. (Such a method makes sense when $\nabla_{x} F$ is much easier to compute than $\nabla g \nabla_{y} F$.) Show that if there exists $\gamma \in(0,1)$ such that

$$
\left\|\nabla g(x) \nabla_{y} F(x, g(x))\right\| \leq \gamma\left\|\nabla_{x} F(x, g(x))\right\|, \quad \forall x \in \Re^{n}
$$

then the method is convergent in the sense that all limit points of the sequences that it generates are stationary points of $f$.
(b) Consider the constrained minimization problem

$$
\begin{aligned}
& \operatorname{minimize} f(x, y) \\
& \text { subject to } h(x, y)=0
\end{aligned}
$$

where $f: \Re^{n+m} \rightarrow \Re$ and $h: \Re^{n+m} \rightarrow \Re^{m}$ are continuously differentiable functions of the two arguments $x \in \Re^{n}$ and $y \in \Re^{m}$. Consider also a method of the form

$$
x^{k+1}=x^{k}-\alpha^{k} \nabla_{x} f\left(x^{k}, y^{k}\right)
$$

where $y^{k}$ is a solution of $h\left(x^{k}, y\right)=0$, viewed as a system of $m$ equations in the unknown vector $y$, and $\alpha^{k}$ is chosen by the minimization rule or the Armijo rule. Formulate conditions that guarantee that this method is convergent.

# 1.3 GRADIENT METHODS - RATE OF CONVERGENCE 

The second major issue regarding gradient methods relates to the rate (or speed) of convergence of the generated sequences $\left\{x^{k}\right\}$. The mere fact that $\left\{x^{k}\right\}$ converges to a stationary point $x^{*}$ will be of little practical value unless the points $x^{k}$ are reasonably close to $x^{*}$ after relatively few iterations. Thus, the study of the rate of convergence provides what are often the dominant criteria for selecting one algorithm in favor of others for solving a particular problem.

## Approaches for Rate of Convergence Analysis

There are several approaches towards quantifying the rate of convergence of nonlinear programming algorithms. We will discuss briefly three possibilities and then concentrate on the third.
(a) Computational complexity approach: Here we try to estimate the number of elementary operations needed by a given method to find

---

an optimal solution exactly or within an $\epsilon$-tolerance. Usually, this approach provides worst-case estimates, that is, upper bounds on the number of required operations over a class of problems of given dimension and type (e.g. linear, convex, etc.). These estimates may also involve parameters such as the distance of the starting point from the optimal solution set, etc.
(b) Informational complexity approach: One difficulty with the computational complexity approach is that for a diverse class of problems, it is often difficult or meaningless to quantify the amount of computation needed for a single function or gradient evaluation. For example, in estimating the computational complexity of the gradient method applied to the entire class of differentiable convex functions, how are we to compare the overhead for finding the stepsize and for updating the $x$ vector with the work needed to compute the cost function value and its gradient? The informational complexity approach, which is discussed in detail in [NeY83] and [TrW80], bypasses this difficulty by estimating the number of function (and possibly gradient) evaluations needed to find an exact or approximately optimal solution (as opposed to the number of necessary computational operations). In other respects, the informational and computational complexity approaches are similar.
(c) Local analysis: In this approach we focus on the local behavior of the method in a neighborhood of an optimal solution. Local analysis can describe quite accurately the behavior of a method near the solution by using Taylor series approximations, but ignores entirely the behavior of the method when far from the solution.

The main potential advantage of the computational and informational complexity approaches is that they provide information about the method's progress when far from the eventual limit. Unfortunately, however, this information is usually pessimistic as it accounts for the worst possible problem instance within the class considered. This has resulted in some striking discrepancies between the theoretical model predictions and practical realworld observations. For example, the most widely used linear programming method, the simplex method, is categorized as a "bad" method by worstcase complexity analysis, because it performs very poorly on some specially constructed examples, which, however, are highly unlikely in practice. On the other hand, the ellipsoid method of Khachiyan [Kha79] (see [BGT81] for a survey), which was the first linear programming method with a polynomial complexity bound, is categorized as much better than the simplex method by worst-case complexity analysis, even though it performs very poorly on most practical linear programs.

The computational complexity approach has received considerable attention in the context of interior point methods. These methods, discussed

---

in Sections 2.6, 4.2, and 4.4, were primarily motivated by Karmarkar's development of a linear programming algorithm with a polynomial complexity bound that was more favorable than the one of the ellipsoid method [Kar84]. It turned out, however, that the worst-case predictions for the required number of iterations of these methods were off by many orders of magnitude from the practically observed number of iterations. Furthermore, the interior point methods that perform best in practice have poor worst-case complexity, while the ones with the best complexity bounds are very slow in practice.

The local analysis approach, which will be adopted exclusively in this text, has enjoyed considerable success in predicting the behavior of various methods near nonsingular local minima where the cost function can be well approximated by a quadratic. However, the local analysis approach also has some important drawbacks, the most important of which is that it does not account for the rate of progress in the initial iterations. Nonetheless, in many practical situations this is not a serious omission because progress is fast in the initial iterations and slows down only in the limit (the reasons for this seem hard to understand; they are problemdependent). Furthermore, often in practice, starting points that are near a solution are easily obtainable by a combination of heuristics and experience. in which case local analysis becomes more meaningful.

Local analysis is not very helpful for problems which either involve singular local minima or which are difficult in the sense that the principal methods take many iterations to get near their solution where local analysis applies. It may be said that at present there is little theory and experience to help a practitioner who is faced with such a problem.

# 1.3.1 The Local Analysis Approach 

We now formalize the basic ingredients of our local rate of convergence analysis approach. These are:
(a) We restrict attention to sequences $\left\{x^{k}\right\}$ that converge to a unique limit point $x^{*}$.
(b) Rate of convergence is evaluated in terms of an error function $e$ : $\Re^{n} \mapsto \Re$ satisfying $e(x) \geq 0$ for all $x \in \Re^{n}$ and $e\left(x^{*}\right)=0$. Typical choices are the Euclidean distance

$$
e(x)=\left\|x-x^{*}\right\|
$$

and the cost difference

$$
e(x)=\left|f(x)-f\left(x^{*}\right)\right|
$$

(c) Our analysis is asymptotic, that is, we look at the rate of convergence of the tail of the error sequence $\left\{e\left(x^{k}\right)\right\}$.

---

(d) The generated error sequence $\left\{e\left(x^{k}\right)\right\}$ is compared with some "standard" sequences. In our case, we compare $\left\{e\left(x^{k}\right)\right\}$ with the geometric progression

$$
\beta^{k}, \quad k=0,1, \ldots
$$

where $\beta \in(0,1)$ is some scalar. In particular, we say that $\left\{e\left(x^{k}\right)\right\}$ converges linearly or geometrically, if there exist $q>0$ and $\beta \in(0,1)$ such that for all $k$

$$
e\left(x^{k}\right) \leq q \beta^{k}
$$

It is possible to show that linear convergence is obtained if for some $\beta \in(0,1)$ we have

$$
\limsup _{k \rightarrow \infty} \frac{e\left(x^{k+1}\right)}{e\left(x^{k}\right)} \leq \beta
$$

that is, asymptotically, the error is dropping by a factor of at least $\beta$ at each iteration (see Exercise 3.6, which gives several additional convergence rate characterizations). If for every $\beta \in(0,1)$, there exists $q$ such that the condition $e\left(x^{k}\right) \leq q \beta^{k}$ holds for all $k$, we say that $\left\{e\left(x^{k}\right)\right\}$ converges superlinearly. This is true in particular, if

$$
\limsup _{k \rightarrow \infty} \frac{e\left(x^{k+1}\right)}{e\left(x^{k}\right)}=0
$$

To quantify further the notion of superlinear convergence, we may compare $\left\{e\left(x^{k}\right)\right\}$ with the sequence

$$
(\beta)^{p^{k}}, \quad k=0,1, \ldots
$$

where $\beta \in(0,1)$, and $p>1$ are some scalars. This sequence converges much faster than a geometric progression. We say that $\left\{e\left(x^{k}\right)\right\}$ converges at least superlinearly with order $p$, if there exist $q>0$, $\beta \in(0,1)$, and $p>1$ such that for all $k$

$$
e\left(x^{k}\right) \leq q(\beta)^{p^{k}}
$$

The case where $p=2$ is referred to as quadratic convergence. It is possible to show that superlinear convergence with order $p$ is obtained if

$$
\limsup _{k \rightarrow \infty} \frac{e\left(x^{k+1}\right)}{e\left(x^{k}\right)^{p}}<\infty
$$

or equivalently, $e\left(x^{k+1}\right)=O\left(e\left(x^{k}\right)^{p}\right)$; see Exercise 3.7.
Most optimization algorithms that are of interest in practice produce sequences converging either linearly or superlinearly, at least when they converge to nonsingular local minima. Linear convergence is a fairly satisfactory rate of convergence for nonlinear programming algorithms, provided the factor $\beta$ of the associated geometric progression is not too close

---

to unity. Several nonlinear programming algorithms converge superlinearly for particular classes of problems. Newton's method is an important example, as we will see in the present section and also in Section 1.4. For convergence to singular local minima, slower than linear convergence rate is quite common.

# 1.3.2 The Role of the Condition Number 

Many of the important convergence rate characteristics of gradient methods reveal themselves when the cost function is quadratic. To see why, assume that a gradient method is applied to minimization of a twice continuously differentiable function function $f: \Re^{n} \mapsto \Re$, and it generates a sequence $\left\{x^{k}\right\}$ converging to a nonsingular local minimum $x^{*}$. By Taylor's theorem we have

$$
f(x)=f\left(x^{*}\right)+\frac{1}{2}\left(x-x^{*}\right)^{\prime} \nabla^{2} f\left(x^{*}\right)\left(x-x^{*}\right)+o\left(\left\|x-x^{*}\right\|^{2}\right)
$$

Therefore, since $\nabla^{2} f\left(x^{*}\right)$ is positive definite, $f$ can be accurately approximated near $x^{*}$ by the quadratic function

$$
f\left(x^{*}\right)+\frac{1}{2}\left(x-x^{*}\right)^{\prime} \nabla^{2} f\left(x^{*}\right)\left(x-x^{*}\right)
$$

We thus expect that asymptotic convergence rate results obtained for the quadratic cost case have direct analogs for the general case. This conjecture can indeed be established by rigorous analysis and has been substantiated by extensive numerical experimentation. For this reason, we take the positive definite quadratic case as our point of departure. We subsequently discuss what happens when $\nabla^{2} f\left(x^{*}\right)$ is not positive definite, in which case an analysis based on a quadratic model is inadequate.

## Convergence Rate of Steepest Descent for Quadratic Functions

Suppose that the cost function $f$ is quadratic with positive definite Hessian $Q$. We may assume without loss of generality that $f$ is minimized at $x^{*}=0$ and that $f\left(x^{*}\right)=0$ [otherwise we can use the change of variables $y=x-x^{*}$ and subtract the constant $f\left(x^{*}\right)$ from $\left.f(x)\right]$. Thus we have

$$
f(x)=\frac{1}{2} x^{\prime} Q x, \quad \nabla f(x)=Q x, \quad \nabla^{2} f(x)=Q
$$

The steepest descent method takes the form

$$
x^{k+1}=x^{k}-\alpha^{k} \nabla f\left(x^{k}\right)=\left(I-\alpha^{k} Q\right) x^{k}
$$

Therefore, we have

$$
\left\|x^{k+1}\right\|^{2}=x^{k^{\prime}}\left(I-\alpha^{k} Q\right)^{2} x^{k}
$$

---

Since by Prop. A.18(b) of Appendix A, we have for all $x \in \Re^{n}$

$$
x^{\prime}\left(I-\alpha^{k} Q\right)^{2} x \leq\left(\text { maximum eigenvalue of }\left(I-\alpha^{k} Q\right)^{2}\right)\|x\|^{2}
$$

we obtain

$$
\left\|x^{k+1}\right\|^{2} \leq\left(\text { maximum eigenvalue of }\left(I-\alpha^{k} Q\right)^{2}\right)\left\|x^{k}\right\|^{2}
$$

Using Prop. A. 13 of Appendix A, it can be seen that the eigenvalues of $\left(I-\alpha^{k} Q\right)^{2}$ are equal to $\left(1-\alpha^{k} \lambda_{i}\right)^{2}$, where $\lambda_{i}$ are the eigenvalues of $Q$. Therefore, we have

$$
\text { maximum eigenvalue of }\left(I-\alpha^{k} Q\right)^{2}=\max \left\{\left(1-\alpha^{k} m\right)^{2},\left(1-\alpha^{k} M\right)^{2}\right\}
$$

where

$$
\begin{aligned}
m: & \text { smallest eigenvalue of } Q \\
M: & \text { largest eigenvalue of } Q
\end{aligned}
$$

It follows that for $x^{k} \neq 0$, we have

$$
\frac{\left\|x^{k+1}\right\|}{\left\|x^{k}\right\|} \leq \max \left\{\left|1-\alpha^{k} m\right|,\left|1-\alpha^{k} M\right|\right\}
$$

It can be seen that if $\left|1-\alpha^{k} m\right| \geq\left|1-\alpha^{k} M\right|$, this inequality holds as an equation if $x^{k}$ is proportional to an eigenvector corresponding to $m$. Otherwise, if $\left|1-\alpha^{k} m\right|<\left|1-\alpha^{k} M\right|$, the inequality holds as an equation if $x^{k}$ is proportional to an eigenvector corresponding to $M$.

Figure 1.3.1 illustrates the convergence rate bound of Eq. (3.1) as a function of the stepsize $\alpha^{k}$. It can be seen that the value of $\alpha^{k}$ that minimizes the bound is

$$
\alpha^{*}=\frac{2}{M+m}
$$

in which case

$$
\frac{\left\|x^{k+1}\right\|}{\left\|x^{k}\right\|} \leq \frac{M-m}{M+m}
$$

This is the best convergence rate bound for steepest descent with constant stepsize.

There is another interesting convergence rate result, which holds when $\alpha^{k}$ is chosen by the line minimization rule. This result quantifies the rate at which the cost decreases and has the form

$$
\frac{f\left(x^{k+1}\right)}{f\left(x^{k}\right)} \leq\left(\frac{M-m}{M+m}\right)^{2}
$$

The above inequality is verified in Prop. 1.3.1, given in the next subsection, where we collect and prove the more formal results of this section. It can

---

![[chapter_21_p60_img19.jpeg]]

Figure 1.3.1. Illustration of the convergence rate bound $\left\|x^{k+1}\right\| /\left\|x^{k}\right\| \leq \max \{|1-$ $\alpha m|,|1-\alpha M|\}$ for steepest descent. The bound is minimized for $\alpha$ such that $1-\alpha m=\alpha M-1$, that is, for $\alpha=2 /(M+m)$.
be shown that the inequality is sharp in the sense that given any $Q$, there is a starting point $x^{0}$ such that this inequality holds as an equation for all $k$ (see Fig. 1.3.2).

The ratio $M / m$ is called the condition number of $Q$, and problems where $M / m$ is large are referred as ill-conditioned. Such problems are characterized by very elongated elliptical level sets. The steepest descent method converges slowly for these problems as indicated by the convergence rate bounds of Eqs. (3.1) and (3.2), and as illustrated in Fig. 1.3.2.

# Scaling and Steepest Descent 

Consider now the more general method

$$
x^{k+1}=x^{k}-\alpha^{k} D^{k} \nabla f\left(x^{k}\right)
$$

where $D^{k}$ is positive definite and symmetric: most of the gradient methods of interest have this form as discussed in Section 1.2. It turns out that we may view this iteration as a scaled version of steepest descent. In particular, this iteration is just steepest descent applied in a different coordinate system, which depends on $D^{k}$.

Indeed, let

$$
S=\left(D^{k}\right)^{1 / 2}
$$

denote the positive definite square root of $D^{k}$ (cf. Prop. A. 21 in Appendix A), and consider a transformation of variables defined by

$$
x=S y
$$

---

![[chapter_21_p61_img20.jpeg]]

Figure 1.3.2. Example showing that the convergence rate bound

$$
\frac{f\left(x^{k+1}\right)}{f\left(x^{k}\right)} \leq\left(\frac{M-m}{M+m}\right)^{2}
$$

is sharp for the steepest descent method with the line minimization rule. Consider the quadratic function

$$
f(x)=\frac{1}{2} \sum_{i=1}^{n} \lambda_{i} x_{i}^{2}
$$

where $0<m=\lambda_{1} \leq \lambda_{2} \leq \cdots \leq \lambda_{n}=M$. Any positive definite quadratic function can be put into this form by transformation of variables. Consider the starting point

$$
x^{0}=\left(m^{-1}, 0, \ldots, 0, M^{-1}\right)^{\prime}
$$

and apply the steepest descent method $x^{k+1}=x^{k}-\alpha^{k} \nabla f\left(x^{k}\right)$ with $\alpha^{k}$ chosen by the line minimization rule. We have $\nabla f\left(x^{0}\right)=(1,0, \ldots, 0,1)^{\prime}$ and it can be verified that the minimizing stepsize is $\alpha^{0}=2 /(M+m)$. Thus we obtain $x_{1}^{1}=1 / m-2 /(M+m), x_{n}^{1}=1 / M-2 /(M+m), x_{i}^{1}=0$ for $i=2, \ldots, n-1$. Therefore,

$$
x^{1}=\left(\frac{M-m}{M+m}\right)\left(m^{-1}, 0, \ldots, 0,-M^{-1}\right)^{\prime}
$$

and, we can verify by induction that for all $k$,

$$
x^{2 k}=\left(\frac{M-m}{M+m}\right)^{2 k} x^{0}, \quad x^{2 k+1}=\left(\frac{M-m}{M+m}\right)^{2 k} x^{1}
$$

Thus, there exist starting points on the plane of points $x$ of the form $x=$ $\left(\xi_{1}, 0, \ldots, 0, \xi_{n}\right)^{\prime}, \xi_{1} \in \Re, \xi_{n} \in \Re$, in fact two lines shown in the figure, for which steepest descent converges in a way that the inequality

$$
\frac{f\left(x^{k+1}\right)}{f\left(x^{k}\right)} \leq\left(\frac{M-m}{M+m}\right)^{2}
$$

is satisfied as an equation at each iteration.

---

Then, in the space of $y$, the problem is written as

$$
\begin{aligned}
& \operatorname{minimize} h(y) \equiv f(S y) \\
& \text { subject to } y \in \Re^{n}
\end{aligned}
$$

The steepest descent method for this problem takes the form

$$
y^{k+1}=y^{k}-\alpha^{k} \nabla h\left(y^{k}\right)
$$

Multiplying with $S$, we obtain

$$
S y^{k+1}=S y^{k}-\alpha^{k} S \nabla h\left(y^{k}\right)
$$

By passing back to the space of $x$, using the relations

$$
x^{k}=S y^{k}, \quad \nabla h\left(y^{k}\right)=S \nabla f\left(x^{k}\right), \quad S^{2}=D^{k}
$$

we obtain

$$
x^{k+1}=x^{k}-\alpha^{k} D^{k} \nabla f\left(x^{k}\right)
$$

Thus the above gradient iteration is nothing but the steepest descent method (3.5) in the space of $y$.

We now apply the convergence rate results for steepest descent to the scaled iteration $y^{k+1}=y^{k}-\alpha^{k} \nabla h\left(y^{k}\right)$, obtaining

$$
\frac{\left\|y^{k+1}\right\|}{\left\|y^{k}\right\|} \leq \max \left\{\left|1-\alpha^{k} m^{k}\right|,\left|1-\alpha^{k} M^{k}\right|\right\}
$$

and

$$
\frac{h\left(y^{k+1}\right)}{h\left(y^{k}\right)} \leq\left(\frac{M^{k}-m^{k}}{M^{k}+m^{k}}\right)^{2}
$$

[cf. the convergence rate bounds (3.1) and (3.2), respectively], where $m^{k}$ and $M^{k}$ are the smallest and largest eigenvalues of the Hessian $\nabla^{2} h(y)$, which is equal to $S \nabla^{2} f(x) S=\left(D^{k}\right)^{1 / 2} Q\left(D^{k}\right)^{1 / 2}$. Using the equation $y^{k}=$ $\left(D^{k}\right)^{-1 / 2} x^{k}$ to pass back to the space of $x$, we obtain the convergence rate bounds

$$
\frac{x^{k+1^{\prime}}\left(D^{k}\right)^{-1} x^{k+1}}{x^{k^{\prime}}\left(D^{k}\right)^{-1} x^{k}} \leq \max \left\{\left(1-\alpha^{k} m^{k}\right)^{2},\left(1-\alpha^{k} M^{k}\right)^{2}\right\}
$$

and

$$
\frac{f\left(x^{k+1}\right)}{f\left(x^{k}\right)} \leq\left(\frac{M^{k}-m^{k}}{M^{k}+m^{k}}\right)^{2}
$$

where

$$
m^{k}: \text { smallest eigenvalue of }\left(D^{k}\right)^{1 / 2} Q\left(D^{k}\right)^{1 / 2}
$$

---

$M^{k}:$ largest eigenvalue of $\left(D^{k}\right)^{1 / 2} Q\left(D^{k}\right)^{1 / 2}$.
The stepsize that minimizes the right-hand side bound of Eq. (3.7) is

$$
\frac{2}{M^{k}+m^{k}}
$$

The important point is that if $M^{k} / m^{k}$ is much larger than unity, the convergence rate can be very slow, even if an optimal stepsize is used. Furthermore, we see that it is desirable to choose $D^{k}$ as close as possible to $Q^{-1}$, so that $\left(D^{k}\right)^{1 / 2}$ is close to $Q^{-1 / 2}$ (cf. Prop. A. 21 in Appendix A) and $M^{k} \approx m^{k} \approx 1$. Note that if $D^{k}$ is so chosen, Eq. (3.9) shows that the stepsize $\alpha=1$ is near optimal.

# Diagonal Scaling 

Many practical problems are ill-conditioned because of poor relative scaling of the optimization variables. By this we mean that the units in which the variables are expressed are incongruent in the sense that single unit changes of different variables have disproportionate effects on the cost.

As an example, consider a financial problem with two variables, investment denoted $x_{1}$ and expressed in dollars, and interest rate denoted $x_{2}$ and expressed in percentage points. If the effect on the cost function $f$ due to a million dollar increment of investment is comparable to the effect due to a percentage point increment of interest rate, then the condition number will be of the order of $10^{12}$ !! [This rough calculation is based on estimating the condition number by the ratio

$$
\frac{\partial^{2} f\left(x_{1}, x_{2}\right)}{\left(\partial x_{2}^{2}\right)} / \frac{\partial^{2} f\left(x_{1}, x_{2}\right)}{\left(\partial x_{1}\right)^{2}}
$$

approximating the second partial derivatives by the finite difference formulas

$$
\begin{aligned}
& \frac{\partial^{2} f\left(x_{1}, x_{2}\right)}{\left(\partial x_{1}\right)^{2}} \approx \frac{f\left(x_{1}+h_{1}, x_{2}\right)+f\left(x_{1}-h_{1}, x_{2}\right)-2 f\left(x_{1}, x_{2}\right)}{h_{1}^{2}} \\
& \frac{\partial^{2} f\left(x_{1}, x_{2}\right)}{\left(\partial x_{2}\right)^{2}} \approx \frac{f\left(x_{1}, x_{2}+h_{2}\right)+f\left(x_{1}, x_{2}-h_{2}\right)-2 f\left(x_{1}, x_{2}\right)}{h_{2}^{2}}
\end{aligned}
$$

and using the relations $f\left(x_{1}+h_{1}, x_{2}\right) \approx f\left(x_{1}, x_{2}+h_{2}\right), f\left(x_{1}-h_{1}, x_{2}\right) \approx$ $f\left(x_{1}, x_{2}-h_{2}\right)$, and $h_{1}=10^{6}, h_{2}=1$, which express the comparability of the effects of a million dollar investment increment and an interest rate percentage point increment.]

The ill-conditioning in such problems can be significantly alleviated by changing the units in which the optimization variables are expressed, which

---

amounts to diagonal scaling of the variables. By this, we mean working in a new coordinate system of a vector $y$ related to $x$ by a transformation,

$$
x=S y
$$

where $S$ is a diagonal matrix. In the absence of further information, a reasonable choice of $S$ is one that makes all the diagonal elements of the Hessian of the cost

$$
S \nabla^{2} f(x) S
$$

in the $y$-coordinate system approximately equal to unity. For this, we must have

$$
s_{i} \approx\left(\frac{\partial^{2} f(x)}{\left(\partial x_{i}\right)^{2}}\right)^{-1 / 2}
$$

where $s_{i}$ is the $i$ th diagonal element of $S$. As discussed earlier, we may express any gradient algorithm in the space of variables $y$ as a gradient algorithm in the space of variables $x$. In particular, steepest descent in the $y$-coordinate system, when translated in the $x$-coordinate system, yields the diagonally scaled steepest descent method

$$
x^{k+1}=x^{k}-\alpha^{k} D^{k} \nabla f\left(x^{k}\right)
$$

where

$$
D^{k}=\left(\begin{array}{cccccccc}
d_{1}^{k} & 0 & 0 & \cdots & 0 & 0 & 0 \\
0 & d_{2}^{k} & 0 & \cdots & 0 & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \\
0 & 0 & 0 & \cdots & 0 & d_{n-1}^{k} & 0 \\
0 & 0 & 0 & \cdots & 0 & 0 & d_{n}^{k}
\end{array}\right)
$$

and

$$
d_{i}^{k} \approx\left(\frac{\partial^{2} f\left(x^{k}\right)}{\left(\partial x_{i}\right)^{2}}\right)^{-1}
$$

This method is also valid for nonquadratic problems as long as $d_{i}^{k}$ are chosen to be positive. It is not guaranteed to improve the convergence rate of steepest descent, but it is simple and often surprisingly effective in practice.

# Nonquadratic Cost Functions 

It is possible to show that our main conclusions on rate of convergence carry over to the nonquadratic case for sequences converging to nonsingular local minima.

Let $f$ be twice continuously differentiable and consider the gradient method

$$
x^{k+1}=x^{k}-\alpha^{k} D^{k} \nabla f\left(x^{k}\right)
$$

---

where $D^{k}$ is positive definite and symmetric. Consider a generated sequence $\left\{x^{k}\right\}$, and assume that

$$
x^{k} \rightarrow x^{*}, \quad \nabla f\left(x^{*}\right)=0, \quad \nabla^{2} f\left(x^{*}\right): \text { positive definite }
$$

and that $x^{k} \neq x^{*}$ for all $k$. Then, denoting

$$
\begin{aligned}
& m^{k}: \text { smallest eigenvalue of }\left(D^{k}\right)^{1 / 2} \nabla^{2} f\left(x^{k}\right)\left(D^{k}\right)^{1 / 2} \\
& M^{k}: \text { largest eigenvalue of }\left(D^{k}\right)^{1 / 2} \nabla^{2} f\left(x^{k}\right)\left(D^{k}\right)^{1 / 2}
\end{aligned}
$$

it is possible to show the following:
(a) There holds

$$
\begin{aligned}
\limsup _{k \rightarrow \infty} & \frac{\left(x^{k+1}-x^{*}\right)^{\prime}\left(D^{k}\right)^{-1}\left(x^{k+1}-x^{*}\right)}{\left(x^{k}-x^{*}\right)^{\prime}\left(D^{k}\right)^{-1}\left(x^{k}-x^{*}\right)} \\
& =\limsup _{k \rightarrow \infty} \max \left\{\left|1-\alpha^{k} m^{k}\right|^{2},\left|1-\alpha^{k} M^{k}\right|^{2}\right\}
\end{aligned}
$$

(b) If $\alpha^{k}$ is chosen by the line minimization rule, there holds

$$
\limsup _{k \rightarrow \infty} \frac{f\left(x^{k+1}\right)-f\left(x^{*}\right)}{f\left(x^{k}\right)-f\left(x^{*}\right)} \leq \limsup _{k \rightarrow \infty}\left(\frac{M^{k}-m^{k}}{M^{k}+m^{k}}\right)^{2}
$$

The proof of these facts essentially involves a repetition of the proofs for the quadratic case. However, the details are complicated and tedious and will not be given.

From Eq. (3.12), we see that if $D^{k}$ converges to some positive definite matrix as $x^{k} \rightarrow x^{*}$, the sequence $\left\{f\left(x^{k}\right)\right\}$ converges to $f\left(x^{*}\right)$ linearly. When

$$
D^{k} \rightarrow \nabla^{2} f\left(x^{*}\right)^{-1}
$$

we have $\lim _{k \rightarrow \infty} M^{k}=\lim _{k \rightarrow \infty} m^{k}=1$ and Eq. (3.12) shows that the convergence rate of $\left\{f\left(x^{k}\right)\right\}$ is superlinear. A somewhat more general version of this result for the case of the Armijo rule is given in the Prop. 1.3.2, which is given in the next subsection. In particular, it is shown that if the direction

$$
d^{k}=-D^{k} \nabla f\left(x^{k}\right)
$$

approaches asymptotically the Newton direction $-\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1} \nabla f\left(x^{k}\right)$ and the Armijo rule is used with initial stepsize equal to one, the rate of convergence is superlinear.

There is a consistent theme that emerges from our analysis, namely that to achieve asymptotically fast convergence of the gradient method

$$
x^{k+1}=x^{k}-\alpha^{k} D^{k} \nabla f\left(x^{k}\right)
$$

---

one should try to choose the matrices $D^{k}$ as close as possible to $\left(\nabla^{2} f\left(x^{*}\right)\right)^{-1}$ so that the maximum and minimum eigenvalues of $\left(D^{k}\right)^{1 / 2} \nabla^{2} f\left(x^{*}\right)\left(D^{k}\right)^{1 / 2}$ satisfy $M^{k} \approx 1$ and $m^{k} \approx 1$. Furthermore, when $D^{k}$ is so chosen, the initial stepsize $s=1$ is a good choice for the Armijo rule and other related rules, or as a starting point for one-dimensional minimization procedures used in minimization stepsize rules. This finding has been supported by extensive numerical experience and is one of the most reliable guidelines for selecting and designing optimization algorithms for unconstrained problems. Note, however, that this guideline is valid only for problems where the cost function is twice differentiable and has positive definite Hessian near the points of interest. We discuss next problems where this condition is not satisfied.

# Singular and Difficult Problems 

We now consider problems where the Hessian matrix either does not exist or is not positive definite at or near local minima of interest. Expressed mathematically, there are local minima $x^{*}$ and directions $d$ such that the slope of $f$ along $d$, which is $\nabla f\left(x^{*}+\alpha d\right)^{\prime} d$, changes very slowly or very rapidly with $\alpha$, that is, either

$$
\lim _{\alpha \rightarrow 0} \frac{\nabla f\left(x^{*}+\alpha d\right)^{\prime} d-\nabla f\left(x^{*}\right)^{\prime} d}{\alpha}=0
$$

or

$$
\lim _{\alpha \rightarrow 0} \frac{\nabla f\left(x^{*}+\alpha d\right)^{\prime} d-\nabla f\left(x^{*}\right)^{\prime} d}{\alpha}=\infty
$$

The case of Eq. (3.13) is characterized by flatness of the cost along the direction $d$; large excursions from $x^{*}$ along $d$ produce small changes in cost. In the case of Eq. (3.14) the reverse is true; the cost rises steeply along $d$. An example is the function

$$
f\left(x_{1}, x_{2}\right)=\left|x_{1}\right|^{4}+\left|x_{2}\right|^{3 / 2}
$$

where for the minimum $x^{*}=(0,0)$, Eq. (3.13) holds along the direction $d=(1,0)$ and Eq. (3.14) holds along the direction $d=(0,1)$. Gradient methods that use directions that are comparable in size to the gradient may require very large stepsizes in the case of Eq. (3.13) and very small stepsizes in the case of Eq. (3.14). This suggests potential difficulties in the implementation of a good stepsize rule; certainly a constant stepsize does not look like an attractive possibility. Furthermore, in the Armijo rule, the initial stepsize should not be taken constant; it should be adjusted according to a suitable scheme, although designing such a scheme may not be easy.

From the point of view of speed of convergence one may view the cases of Eqs. (3.13) and (3.14) as corresponding to an "infinite condition number," thereby suggesting slower than linear convergence rate for the

---

method of steepest descent. Proposition 1.3.3 of the next subsection quantifies the rate of convergence of gradient methods for the case of a convex function whose gradient satisfies the Lipschitz condition

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|
$$

for some $L$, and all $x$ and $y$ in a neighborhood of $x^{*}$ [this assumption is consistent with the "flat" cost case of Eq. (3.13), but not with the "steep" cost case of Eq. (3.14)]. It is shown in particular that for a gradient method with the minimization rule, we have

$$
f\left(x^{k}\right)-f\left(x^{*}\right)=o(1 / k)
$$

This type of estimate suggests that for many practical singular problems one may be unable to obtain a highly accurate approximation of an optimal solution. In the "steep" cost case where Eq. (3.14) holds for some directions $d$, computational examples suggest that the rate of convergence can be slower than linear for the method of steepest descent, although a formal analysis of this conjecture does not seem to have been published.

It should be noted that problems with singular local minima are not the only ones for which gradient methods may converge slowly. There are problems where a given method may have excellent asymptotic rate of convergence, but its progress when far from the eventual limit can be very slow. A prominent example is when the cost function is continuously differentiable but its Hessian matrix is discontinuous and possibly singular in some regions that are outside a small neighborhood of the solution; such functions arise for example in augmented Lagrangian methods for inequality constrained problems (see Section 4.2). Then the powerful Newton-like methods may require a very large number of iterations to get to the small neighborhood of the eventual limit where their convergence rate is favorable. What happens here is that these methods use second derivative information in sophisticated ways, but this information may be misleading due to the Hessian discontinuities.

Generally, there is a tendency to think that difficult problems should be addressed with sophisticated methods, such as Newton-like methods. This is often true, particularly for problems with nonsingular local minima that are poorly conditioned. However, it is important to realize that often the reverse is true, namely that for problems with "difficult" cost functions and singular local minima, it is best to use simple methods such as (perhaps diagonally scaled) steepest descent with simple stepsize rules such as a constant or a diminishing stepsize. The reason is that methods that use sophisticated descent directions and stepsize rules often rely on assumptions that are likely to be violated in difficult problems. We also note that for difficult problems, it may be helpful to supplement the steepest descent method with features that allow it to deal better with multiple local minima and peculiarities of the cost function. An often useful modification is the heavy ball method, discussed in Exercise 3.9.

---

# 1.3.3 Convergence Rate Results 

We first derive the convergence rate of steepest descent with the minimization stepsize rule when the cost is quadratic.

Proposition 1.3.1: Consider the quadratic function

$$
f(x)=\frac{1}{2} x^{\prime} Q x
$$

where $Q$ is positive definite and symmetric, and the method of steepest descent

$$
x^{k+1}=x^{k}-\alpha^{k} \nabla f\left(x^{k}\right)
$$

where the stepsize $\alpha^{k}$ is chosen according to the minimization rule

$$
f\left(x^{k}-\alpha^{k} \nabla f\left(x^{k}\right)\right)=\min _{\alpha \geq 0} f\left(x^{k}-\alpha \nabla f\left(x^{k}\right)\right)
$$

Then, for all $k$,

$$
f\left(x^{k+1}\right) \leq\left(\frac{M-m}{M+m}\right)^{2} f\left(x^{k}\right)
$$

where $M$ and $m$ are the largest and smallest eigenvalues of $Q$, respectively.

Proof: Let us denote

$$
g^{k}=\nabla f\left(x^{k}\right)=Q x^{k}
$$

The result clearly holds if $g^{k}=0$, so we assume $g^{k} \neq 0$. We first compute the minimizing stepsize $\alpha^{k}$. We have

$$
\frac{d}{d \alpha} f\left(x^{k}-\alpha g^{k}\right)=-g^{k^{\prime}} Q\left(x^{k}-\alpha g^{k}\right)=-g^{k^{\prime}} g^{k}+\alpha g^{k^{\prime}} Q g^{k}
$$

By setting this derivative equal to zero, we obtain

$$
\alpha^{k}=\frac{g^{k^{\prime}} g^{k}}{g^{k^{\prime}} Q g^{k}}
$$

We have, using Eqs. (3.16)-(3.18),

$$
\begin{aligned}
f\left(x^{k+1}\right) & =\frac{1}{2}\left(x^{k}-\alpha^{k} g^{k}\right)^{\prime} Q\left(x^{k}-\alpha^{k} g^{k}\right) \\
& =\frac{1}{2}\left(x^{k^{\prime}} Q x^{k}-2 \alpha^{k} g^{k^{\prime}} Q x^{k}+\left(\alpha^{k}\right)^{2} g^{k^{\prime}} Q g^{k}\right) \\
& =\frac{1}{2}\left(x^{k^{\prime}} Q x^{k}-2 \alpha^{k} g^{k^{\prime}} g^{k}+\left(\alpha^{k}\right)^{2} g^{k^{\prime}} Q g^{k}\right)
\end{aligned}
$$

---

and using Eq. (3.19),

$$
f\left(x^{k+1}\right)=\frac{1}{2}\left(x^{k^{\prime}} Q x^{k}-\frac{\left(g^{k^{\prime}} g^{k}\right)^{2}}{g^{k^{\prime}} Q g^{k}}\right)
$$

Thus, using the fact $f\left(x^{k}\right)=\frac{1}{2} x^{k^{\prime}} Q x^{k}=\frac{1}{2} g^{k^{\prime}} Q^{-1} g^{k}$, we obtain

$$
f\left(x^{k+1}\right)=\left(1-\frac{\left(g^{k^{\prime}} g^{k}\right)^{2}}{\left(g^{k^{\prime}} Q g^{k}\right)\left(g^{k^{\prime}} Q^{-1} g^{k}\right)}\right) f\left(x^{k}\right)
$$

At this point we need the following lemma.

Lemma 3.1: (Kantorovich Inequality) Let $Q$ be a positive definite and symmetric $n \times n$ matrix. Then for any vector $y \in \Re^{n}, y \neq 0$, there holds

$$
\frac{\left(y^{\prime} y\right)^{2}}{\left(y^{\prime} Q y\right)\left(y^{\prime} Q^{-1} y\right)} \geq \frac{4 M m}{(M+m)^{2}}
$$

where $M$ and $m$ are the largest and smallest eigenvalues of $Q$, respectively.

Proof: Let $\lambda_{1}, \ldots, \lambda_{n}$ denote the eigenvalues of $Q$ and assume that

$$
0<m=\lambda_{1} \leq \lambda_{2} \leq \cdots \leq \lambda_{n}=M
$$

Let $S$ be the matrix consisting of the $n$ orthogonal eigenvectors of $Q$, normalized so that they have unit norm (cf. Prop. A. 27 in Appendix A). Then, it can be seen that $S^{\prime} Q S$ is diagonal with diagonal elements $\lambda_{1}, \ldots, \lambda_{n}$. By using if necessary a transformation of the coordinate system that replaces $y$ by $S x$, we may assume that $Q$ is diagonal and that its diagonal elements are $\lambda_{1}, \ldots, \lambda_{n}$. We have for $y=\left(y_{1}, \ldots, y_{n}\right)^{\prime} \neq 0$

$$
\frac{\left(y^{\prime} y\right)^{2}}{\left(y^{\prime} Q y\right)\left(y^{\prime} Q^{-1} y\right)}=\frac{\left(\sum_{i=1}^{n} y_{i}^{2}\right)^{2}}{\left(\sum_{i=1}^{n} \lambda_{i} y_{i}^{2}\right)\left(\sum_{i=1}^{n} \frac{y_{i}^{2}}{\lambda_{i}}\right)}
$$

By letting

$$
\xi_{j}=\frac{y_{j}^{2}}{\sum_{i=1}^{n} y_{i}^{2}}
$$

and by defining

$$
\phi(\xi)=\frac{1}{\sum_{i=1}^{n} \xi_{i} \lambda_{i}}, \quad \psi(\xi)=\sum_{i=1}^{n} \frac{\xi_{i}}{\lambda_{i}}
$$

---

we obtain

$$
\frac{\left(y^{\prime} y\right)^{2}}{\left(y^{\prime} Q y\right)\left(y^{\prime} Q^{-1} y\right)}=\frac{\phi(\xi)}{\psi(\xi)}
$$

Figure 1.3.3 shows that we have

$$
\frac{\phi(\xi)}{\psi(\xi)} \geq \frac{4 \lambda_{1} \lambda_{n}}{\left(\lambda_{1}+\lambda_{n}\right)^{2}}
$$

which proves the desired inequality. Q.E.D.
Returning to the proof of Prop. 1.3.1, we have by using the Kantorovich inequality in Eq. (3.20)

$$
f\left(x^{k+1}\right) \leq\left(1-\frac{4 M m}{(M+m)^{2}}\right) f\left(x^{k}\right)=\left(\frac{M-m}{M+m}\right)^{2} f\left(x^{k}\right)
$$

# Q.E.D. 

The following proposition shows superlinear convergence for methods where $d^{k}$ approaches the Newton direction $-\left(\nabla^{2} f\left(x^{*}\right)\right)^{-1} \nabla f\left(x^{k}\right)$ and the Armijo rule is used.

## Proposition 1.3.2: (Superlinear Convergence of Newton-Like

Methods) Let $f$ be twice continuously differentiable. Consider a sequence $\left\{x^{k}\right\}$ generated by the gradient method $x^{k+1}=x^{k}+\alpha^{k} d^{k}$ and suppose that

$$
x^{k} \rightarrow x^{*}, \quad \nabla f\left(x^{*}\right)=0, \quad \nabla^{2} f\left(x^{*}\right): \text { positive definite. }
$$

Assume further that $\nabla f\left(x^{k}\right) \neq 0$ for all $k$ and

$$
\lim _{k \rightarrow \infty} \frac{\left\|d^{k}+\left(\nabla^{2} f\left(x^{*}\right)\right)^{-1} \nabla f\left(x^{k}\right)\right\|}{\left\|\nabla f\left(x^{k}\right)\right\|}=0
$$

Then, if $\alpha^{k}$ is chosen by means of the Armijo rule with initial stepsize $s=1$ and $\sigma<1 / 2$, we have

$$
\lim _{k \rightarrow \infty} \frac{\left\|x^{k+1}-x^{*}\right\|}{\left\|x^{k}-x^{*}\right\|}=0
$$

Furthermore, there exists an integer $\bar{k} \geq 0$ such that $\alpha^{k}=1$ for all $k \geq \bar{k}$ (i.e., eventually no reduction of the initial stepsize will be taking place).

---

![[chapter_21_p71_img21.jpeg]]

Figure 1.3.3. Proof of the Kantorovich inequality. Consider the function $1 / \lambda$. The scalar $\sum_{i=1}^{n} \xi_{i} \lambda_{i}$ represents, for any $\xi=\left(\xi_{1}, \ldots, \xi_{n}\right)$ with $\xi_{i} \geq 0, \sum_{i=1}^{n} \xi_{i}=$ 1 , a point in the line segment $\left[\lambda_{1}, \lambda_{n}\right]$. Thus, the values $\phi(\xi)=1 / \sum_{i=1}^{n} \xi_{i} \lambda_{i}$ correspond to the thick part of the curve $1 / \lambda$. On the other hand, the value $\psi(\xi)=$ $\sum_{i=1}^{n}\left(\xi_{i} / \lambda_{i}\right)$ is a convex combination of $1 / \lambda_{1}, \ldots, 1 / \lambda_{n}$ and hence corresponds to a point in the shaded area in the figure. For the same vector $\xi$, both $\phi(\xi)$ and $\psi(\xi)$ are represented by points on the same vertical line. Hence,

$$
\frac{\phi(\xi)}{\psi(\xi)} \geq \min _{\lambda_{1} \leq \lambda \leq \lambda_{n}} \frac{\frac{1}{\lambda}}{\frac{\lambda_{1}+\lambda_{n}-\lambda}{\lambda_{1} \lambda_{n}}}
$$

The minimum is attained for $\lambda=\left(\lambda_{1}+\lambda_{n}\right) / 2$ and we obtain

$$
\frac{\phi(\xi)}{\psi(\xi)} \geq \frac{4 \lambda_{1} \lambda_{n}}{\left(\lambda_{1}+\lambda_{n}\right)^{2}}
$$

which is used to show the result.

Proof: We first prove that there exists a $\bar{k} \geq 0$ such that for all $k \geq \bar{k}$,

$$
f\left(x^{k}+d^{k}\right)-f\left(x^{k}\right) \leq \sigma \nabla f\left(x^{k}\right)^{\prime} d^{k}
$$

that is, the unity initial stepsize passes the test of the Armijo rule. By the mean value theorem, we have

$$
f\left(x^{k}+d^{k}\right)-f\left(x^{k}\right)=\nabla f\left(x^{k}\right)^{\prime} d^{k}+\frac{1}{2} d^{k^{\prime}} \nabla^{2} f\left(\bar{x}^{k}\right) d^{k}
$$

where $\bar{x}^{k}$ is a point on the line segment joining $x^{k}$ and $x^{k}+d^{k}$. Thus, it will be sufficient to show that for $k$ sufficiently large, we have

$$
\nabla f\left(x^{k}\right)^{\prime} d^{k}+\frac{1}{2} d^{k} \nabla^{2} f\left(\bar{x}^{k}\right) d^{k} \leq \sigma \nabla f\left(x^{k}\right)^{\prime} d^{k}
$$

---

By defining

$$
p^{k}=\frac{\nabla f\left(x^{k}\right)}{\left\|\nabla f\left(x^{k}\right)\right\|}, \quad q^{k}=\frac{d^{k}}{\left\|\nabla f\left(x^{k}\right)\right\|}
$$

this condition is written

$$
(1-\sigma) p^{k^{\prime}} q^{k}+\frac{1}{2} q^{k^{\prime}} \nabla^{2} f\left(\bar{x}^{k}\right) q^{k} \leq 0
$$

From Eq. (3.22), we have $q^{k}-\left(\nabla^{2} f\left(x^{*}\right)\right)^{-1} p^{k} \rightarrow 0$. Since $\nabla^{2} f\left(x^{*}\right)$ is positive definite and $\left\|p^{k}\right\|=1$, it follows that $\left\{q^{k}\right\}$ is a bounded sequence, and in view of $q^{k}=d^{k} /\left\|\nabla f\left(x^{k}\right)\right\|$ and $\nabla f\left(x^{k}\right) \rightarrow 0$, we obtain $d^{k} \rightarrow 0$. Hence, $x^{k}+d^{k} \rightarrow x^{*}$, and it follows that $\bar{x}^{k} \rightarrow x^{*}$ and $\nabla^{2} f\left(\bar{x}^{k}\right) \rightarrow \nabla^{2} f\left(x^{*}\right)$. We now write Eq. (3.22) as

$$
q^{k}=-\left(\nabla^{2} f\left(x^{*}\right)\right)^{-1} p^{k}+\beta^{k}
$$

where $\left\{\beta^{k}\right\}$ denotes a vector sequence with $\beta^{k} \rightarrow 0$. By using the above relation and the fact $\nabla^{2} f\left(\bar{x}^{k}\right) \rightarrow \nabla^{2} f\left(x^{*}\right)$, we may write Eq. (3.24) as

$$
(1-\sigma) p^{k^{\prime}}\left(\nabla^{2} f\left(x^{*}\right)\right)^{-1} p^{k}-\frac{1}{2} p^{k^{\prime}}\left(\nabla^{2} f\left(x^{*}\right)\right)^{-1} p^{k} \geq \gamma^{k}
$$

where $\left\{\gamma^{k}\right\}$ is some scalar sequence with $\gamma^{k} \rightarrow 0$. Thus Eq. (3.24) is equivalent to

$$
\left(\frac{1}{2}-\sigma\right) p^{k^{\prime}}\left(\nabla^{2} f\left(x^{*}\right)\right)^{-1} p^{k} \geq \gamma^{k}
$$

Since $1 / 2>\sigma,\left\|p^{k}\right\|=1$, and $\nabla^{2} f\left(x^{*}\right)$ is positive definite, the above relation holds for sufficiently large $k$. Thus, the unity initial stepsize is acceptable for sufficiently large $k$, as desired.

To complete the proof, we note that from Eq. (3.22), we have

$$
d^{k}+\left(\nabla^{2} f\left(x^{*}\right)\right)^{-1} \nabla f\left(x^{k}\right)=\left\|\nabla f\left(x^{k}\right)\right\| \delta^{k}
$$

where $\delta^{k}$ is some vector sequence with $\delta^{k} \rightarrow 0$. From Taylor's theorem we obtain

$$
\nabla f\left(x^{k}\right)=\nabla^{2} f\left(x^{*}\right)\left(x^{k}-x^{*}\right)+o\left(\left\|x^{k}-x^{*}\right\|\right)
$$

from which

$$
\begin{aligned}
\left(\nabla^{2} f\left(x^{*}\right)\right)^{-1} \nabla f\left(x^{k}\right) & =x^{k}-x^{*}+o\left(\left\|x^{k}-x^{*}\right\|\right) \\
\left\|\nabla f\left(x^{k}\right)\right\| & =O\left(\left\|x^{k}-x^{*}\right\|\right)
\end{aligned}
$$

Using the above two relations in Eq. (3.25), we obtain

$$
d^{k}+x^{k}-x^{*}=o\left(\left\|x^{k}-x^{*}\right\|\right)
$$

Since for sufficiently large $k$ we have $d^{k}+x^{k}=x^{k+1}$, Eq. (3.26) yields

$$
x^{k+1}-x^{*}=o\left(\left\|x^{k}-x^{*}\right\|\right)
$$

---

from which

$$
\lim _{k \rightarrow \infty} \frac{\left\|x^{k+1}-x^{*}\right\|}{\left\|x^{k}-x^{*}\right\|}=\lim _{k \rightarrow \infty} \frac{o\left(\left\|x^{k}-x^{*}\right\|\right)}{\left\|x^{k}-x^{*}\right\|}=0
$$

# Q.E.D. 

Note that the equation $\lim _{k \rightarrow \infty}\left(\left\|x^{k+1}-x^{*}\right\| /\left\|x^{k}-x^{*}\right\|\right)=0$ implies that $\left\{\left\|x^{k}-x^{*}\right\|\right\}$ converges superlinearly (see Exercise 3.6). In particular, we see that Newton's method, combined with the Armijo rule with unity initial stepsize, converges superlinearly when it converges to a local minimum $x^{*}$ such that $\nabla^{2} f\left(x^{*}\right)$ is positive definite. The capture theorem (Prop. 1.2.5) together with the preceding proposition suggest that Newtonlike methods with the Armijo rule and a unity initial stepsize converge to a local minimum $x^{*}$ such that $\nabla^{2} f\left(x^{*}\right)$ is positive definite, whenever they are started sufficiently close to such a local minimum. The proof of this is left as Exercise 3.2 for the reader.

We finally consider the convergence rate of gradient methods for singular problems whose cost is sufficiently flat for a Lipschitz condition on the gradient to hold.

Proposition 1.3.3: (Convergence Rate of Gradient Methods for Singular Problems) Suppose that the cost function $f$ is convex and its gradient satisfies for some $L$ the Lipschitz condition

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|, \quad \forall x, y \in \Re^{n}
$$

Consider a gradient method $x^{k+1}=x^{k}+\alpha^{k} d^{k}$ where $\alpha^{k}$ is chosen by the minimization rule, and the angle between $d^{k}$ and $\nabla f\left(x^{k}\right)$ is bounded away from 90 degrees, that is, for some $c>0$ and all $k$ we have

$$
\nabla f\left(x^{k}\right)^{\prime} d^{k} \leq-c\left\|\nabla f\left(x^{k}\right)\right\|\left\|d^{k}\right\|
$$

Suppose that the set of global minima $X^{*}$ of $f$ is nonempty and bounded. Then

$$
f\left(x^{k}\right)-f^{*}=o(1 / k)
$$

where $f^{*}=\min _{x} f(x)$ is the optimal value.

Proof: We assume that $\nabla f\left(x^{k}\right) \neq 0$ and therefore also $d^{k} \neq 0$ for all $k$; otherwise the method terminates finitely at a global minimum and the result holds trivially. Let

$$
\tilde{\alpha}^{k}=\frac{\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right|}{L\left\|d^{k}\right\|^{2}}
$$

---

$$
\tilde{x}^{k}=x^{k}+\tilde{\alpha}^{k} d^{k}
$$

By using the descent lemma (Prop. A. 24 in Appendix A), and Eqs. (3.28) and (3.29), we have

$$
\begin{aligned}
f\left(\tilde{x}^{k}\right)-f\left(x^{k}\right) & \leq-\tilde{\alpha}^{k}\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right|+\frac{1}{2}\left(\tilde{\alpha}^{k}\right)^{2} L\left\|d^{k}\right\|^{2} \\
& =\tilde{\alpha}^{k}\left(-\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right|+\frac{1}{2}\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right|\right) \\
& =-\frac{\tilde{\alpha}^{k}}{2}\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right| \\
& =-\frac{\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right|^{2}}{2 L\left\|d^{k}\right\|^{2}} \\
& \leq-\frac{c^{2}\left\|\nabla f\left(x^{k}\right)\right\|^{2}}{2 L}
\end{aligned}
$$

Using this relation together with the fact $f\left(x^{k+1}\right) \leq f\left(\tilde{x}^{k}\right)$, we obtain

$$
f\left(x^{k+1}\right) \leq f\left(x^{k}\right)-\frac{c^{2}\left\|\nabla f\left(x^{k}\right)\right\|^{2}}{2 L}
$$

Since $X^{*}$, the set of global minima of $f$, is nonempty and compact, all the level sets of $f$ are compact (Prop. B. 9 in Appendix B). Thus, $\left\{x^{k}\right\}$ is bounded, and by Prop. 1.2.1, all limit points of $\left\{x^{k}\right\}$ belong to $X^{*}$, and the distance of $x^{k}$ from $X^{*}$, defined by

$$
d\left(x^{k}, X^{*}\right)=\min _{x^{*} \in X^{*}}\left\|x^{k}-x^{*}\right\|
$$

converges to 0 . Using the convexity of $f$, we also have for every global minimum $x^{*}$

$$
f\left(x^{k}\right)-f\left(x^{*}\right) \leq \nabla f\left(x^{k}\right)^{\prime}\left(x^{k}-x^{*}\right) \leq\left\|\nabla f\left(x^{k}\right)\right\| \cdot\left\|x^{k}-x^{*}\right\|
$$

from which, by minimizing over $x^{*} \in X^{*}$,

$$
f\left(x^{k}\right)-f^{*} \leq\left\|\nabla f\left(x^{k}\right)\right\| d\left(x^{k}, X^{*}\right)
$$

Let us denote for all $k$

$$
e^{k}=f\left(x^{k}\right)-f^{*}
$$

Combining Eqs. (3.31) and (3.32), we obtain

$$
e^{k+1} \leq e^{k}-\frac{c^{2}\left(e^{k}\right)^{2}}{2 L d\left(x^{k}, X^{*}\right)^{2}} . \quad \forall k
$$

where we assume without loss of generality that $d\left(x^{k}, X^{*}\right) \neq 0$.
We will show that Eq. (3.33) implies that $e^{k}=o(1 / k)$. Indeed we have

$$
0<e^{k+1} \leq e^{k}\left(1-\frac{c^{2} e^{k}}{2 L d\left(x^{k}, X^{*}\right)^{2}}\right)
$$

---

$$
0<1-\frac{c^{2} e^{k}}{2 L d\left(x^{k}, X^{*}\right)^{2}}
$$

from which

$$
\begin{aligned}
\left(e^{k+1}\right)^{-1} & \geq\left(e^{k}\right)^{-1}\left(1-\frac{c^{2} e^{k}}{2 L d\left(x^{k}, X^{*}\right)^{2}}\right)^{-1} \geq\left(e^{k}\right)^{-1}\left(1+\frac{c^{2} e^{k}}{2 L d\left(x^{k}, X^{*}\right)^{2}}\right) \\
& =\left(e^{k}\right)^{-1}+\frac{c^{2}}{2 L d\left(x^{k}, X^{*}\right)^{2}}
\end{aligned}
$$

Summing this inequality over all $k$, we obtain

$$
e^{k} \leq\left(\left(e^{0}\right)^{-1}+\frac{c^{2}}{2 L} \sum_{i=0}^{k-1} d\left(x^{i}, X^{*}\right)^{-2}\right)^{-1}
$$

or

$$
k e^{k} \leq\left(\frac{1}{k e^{0}}+\frac{c^{2}}{2 L k} \sum_{i=0}^{k-1} d\left(x^{i}, X^{*}\right)^{-2}\right)^{-1}
$$

Since $d\left(x^{i}, X^{*}\right) \rightarrow 0$, we have $d\left(x^{i}, X^{*}\right)^{-2} \rightarrow \infty$ and

$$
\frac{c^{2}}{2 L k} \sum_{i=0}^{k-1} d\left(x^{i}, X^{*}\right)^{-2} \rightarrow \infty
$$

Therefore the right-hand side of Eq. (3.34) tends to 0 , implying that $e^{k}=$ $o(1 / k)$. Q.E.D.

Note that the preceding proof can be modified to cover the case where the Lipschitz condition (3.27) holds within the set $\left\{x \mid f(x) \leq f\left(x^{0}\right)\right\}$. Furthermore, the proof goes through for any stepsize rule for which a relation of the form $f\left(x^{k+1}\right) \leq f\left(x^{k}\right)-\gamma\|\nabla f\left(x^{k}\right)\|^{2}$ can be established for some $\gamma>0$ [cf. Eq. (3.31)]; see Exercise 3.8.

With additional assumptions on the structure of the function $f$ some more precise convergence rate results can be obtained. In particular, if $f$ is convex, has a unique minimum $x^{*}$, and satisfies the following growth condition

$$
f(x)-f\left(x^{*}\right) \geq q\left\|x-x^{*}\right\|^{\beta}, \quad \forall x \text { such that } f(x) \leq f\left(x^{0}\right)
$$

for some scalars $q>0$ and $\beta>2$, it can be shown [Dun81] that for the method of steepest descent with the Armijo rule we have

$$
f\left(x^{k}\right)-f\left(x^{*}\right)=O\left(\frac{1}{k^{\frac{\beta}{\beta-3}}}\right)
$$

---

# E XERCISES 

## 3.1

Estimate the rate of convergence of steepest descent with the line minimization rule when applied to the function of two variables $f(x, y)=x^{2}+1.999 x y+$ $y^{2}$. Find a starting point for which this estimate is sharp (cf. Fig. 1.3.2).

## 3.2

Let $f$ be twice continuously differentiable. Consider a sequence $\left\{x^{k}\right\}$ generated by the gradient method $x^{k+1}=x^{k}+\alpha^{k} d^{k}$ and suppose that $x^{*}$ is a nonsingular local minimum. Assume that, for all $k, \nabla f\left(x^{k}\right) \neq 0$ and $d^{k}=d\left(x^{k}\right)$, where $d(\cdot)$ is a continuous function of $x$ with

$$
\lim _{x \rightarrow x^{*}, \nabla f(x) \neq 0} \frac{\left\|d(x)+\left(\nabla^{2} f(x)\right)^{-1} \nabla f(x)\right\|}{\|\nabla f(x)\|}=0
$$

Furthermore, $\alpha^{k}$ is chosen by means of the Armijo rule with initial stepsize $s=1$ and $\sigma<1 / 2$. Show that there exists an $\epsilon>0$ such that if $\left\|x^{0}-x^{*}\right\|<\epsilon$, then:
(a) $\left\{x^{k}\right\}$ converges to $x^{*}$.
(b) $\alpha^{k}=1$ for all $k$.
(c) $\lim _{k \rightarrow \infty}\left(\left\|x^{k+1}-x^{*}\right\| /\left\|x^{k}-x^{*}\right\|\right)=0$.

## 3.3

Consider a positive definite quadratic problem with Hessian matrix $Q$. Suppose we use scaling with the diagonal matrix whose $i$ th diagonal element is $q_{i i}^{-1}$, where $q_{i i}$ is the $i$ th diagonal element of $Q$. Show that if $Q$ is $2 \times 2$, this diagonal scaling improves the condition number of the problem and the convergence rate of steepest descent. (Note: This need not be true for dimensions higher than 2.)

### 3.4 (Steepest Descent with Errors)

Consider the steepest descent method

$$
x^{k+1}=x^{k}-s\left(\nabla f\left(x^{k}\right)+e^{k}\right)
$$

where $s$ is a constant stepsize, $e^{k}$ is an error satisfying $\left\|e^{k}\right\| \leq \delta$ for all $k$, and $f$ is the positive definite quadratic function

$$
f(x)=\frac{1}{2}\left(x-x^{*}\right)^{\prime} Q\left(x-x^{*}\right)
$$

---

Let

$$
q=\max \{|1-s m|,|1-s M|\}
$$

where
$m$ : smallest eigenvalue of $Q, \quad M$ : largest eigenvalue of $Q$,
and assume that $q<1$. Show that for all $k$, we have

$$
\left\|x^{k}-x^{*}\right\| \leq \frac{s \delta}{1-q}+q^{k}\left\|x^{0}-x^{*}\right\|
$$

# 3.5 

Consider the positive definite quadratic function $f(x)=\frac{1}{2} x^{\prime} Q x$ and the steepest descent method with the stepsize $\alpha^{k}$ chosen by the Goldstein rule. Show that for all $k$,

$$
f\left(x^{k+1}\right) \leq\left(1-\frac{16 \sigma(1-\sigma) M m}{(M+m)^{2}}\right) f\left(x^{k}\right)
$$

Explain why when $\sigma=1 / 2$ this relation yields the result of Prop. 1.3.1. Hint: Use the result of Exercise 2.9.

## $3.6[$ Ber82a]

Consider a scalar sequence $\left\{e^{k}\right\}$ with $e^{k} \geq 0$ for all $k$, and $e^{k} \rightarrow 0$. We say that $\left\{e^{k}\right\}$ converges faster than linearly with convergence ratio $\beta$, where $0<\beta<1$, if for every $\bar{\beta} \in(\beta, 1)$ and $q>0$, there exists $\bar{k}$ such that

$$
e^{k} \leq q \bar{\beta}^{k} . \quad \forall k \geq \bar{k}
$$

We say that $\left\{e^{k}\right\}$ converges slower than linearly with convergence ratio $\beta$, where $0<\beta<1$, if for every $\bar{\beta} \in(\beta, 1)$ and $q>0$, there exists $\bar{k}$ such that

$$
q \bar{\beta}^{k} \leq e^{k}, \quad \forall k \geq \bar{k}
$$

We say that $\left\{e^{k}\right\}$ converges linearly with convergence ratio $\beta$ if it converges both faster and slower than linearly with convergence ratio $\beta$. Show that:
(a) $\left\{e^{k}\right\}$ converges faster than linearly with convergence ratio $\beta$ if and only if

$$
\limsup _{k \rightarrow \infty}\left(e^{k}\right)^{1 / k} \leq \beta
$$

$\left\{e^{k}\right\}$ converges slower than linearly with convergence ratio $\beta$ if and only if

$$
\liminf _{k \rightarrow \infty}\left(e^{k}\right)^{1 / k} \geq \beta
$$

---

$\left\{e^{k}\right\}$ converges linearly with convergence ratio $\beta$ if and only if

$$
\lim _{k \rightarrow \infty}\left(e^{k}\right)^{1 / k}=\beta
$$

(b) Assume that $e^{k} \neq 0$ for all $k$, and denote

$$
\beta_{1}=\liminf _{k \rightarrow \infty} \frac{e^{k+1}}{e^{k}}, \quad \beta_{2}=\limsup _{k \rightarrow \infty} \frac{e^{k+1}}{e^{k}}
$$

Show that if $0<\beta_{1} \leq \beta_{2}<1$, then $\left\{e^{k}\right\}$ converges faster than linearly with convergence ratio $\beta_{2}$ and slower than linearly with convergence ratio $\beta_{1}$. Furthermore, if $\beta_{1}=\beta_{2}=0$, then $\left\{e^{k}\right\}$ converges superlinearly.

# 3.7 

Consider a scalar sequence $\left\{e^{k}\right\}$ with $e^{k}>0$ for all $k$, and $e^{k} \rightarrow 0$. Show that $\left\{e^{k}\right\}$ converges superlinearly with order $p$ if

$$
\limsup _{k \rightarrow \infty} \frac{e^{k+1}}{\left(e^{k}\right)^{p}}<\infty
$$

3.8

Prove the result of Prop. 1.3.3 for the steepest descent case $\left[d^{k}=-\nabla f\left(x^{k}\right)\right]$, and assuming that the stepsize is not chosen by the line minimization rule but is instead a sufficiently small constant.

## 3.9 (The Heavy Ball Method [Pol64])

Consider the following variant of the steepest descent method:

$$
x^{k+1}=x^{k}-\alpha \nabla f\left(x^{k}\right)+\beta\left(x^{k}-x^{k-1}\right), \quad k=1,2, \ldots
$$

where $\alpha$ is a constant positive stepsize and $\beta$ is a scalar with $0<\beta<1$.
(a) Let $f$ be the quadratic function $f(x)=(1 / 2) x^{\prime} Q x+c^{\prime} x$, where $Q$ is positive definite and symmetric, and let $m$ and $M$ be the minimum and the maximum eigenvalues of $Q$, respectively. Show that the method converges linearly to the unique solution if $0<\alpha<2(1+\beta) / M$. Show that with optimal choices of $\alpha$ and $\beta$, the ratio of linear convergence is

$$
\frac{\sqrt{M}-\sqrt{m}}{\sqrt{M}+\sqrt{m}}
$$

---

which, if $m<M$, is faster than the corresponding ratio of the steepest descent method where $\beta=0$ and $\alpha$ is chosen optimally [cf. Eq. (3.1)]. Hint: Consider the iteration

$$
\binom{x^{k+1}}{x^{k}}=\left(\begin{array}{cc}
(1+\beta) I-\alpha Q & -\beta I \\
I & 0
\end{array}\right)\binom{x^{k}}{x^{k-1}}
$$

and show that $v$ is an eigenvalue of the matrix in the above equation if and only if $v+\beta / v$ is equal to $1+\beta-\alpha \lambda$ where $\lambda$ is an eigenvalue of $Q$.
(b) It is generally conjectured that in comparison to steepest descent, the method is less prone to getting trapped at "shallow" local minima, and tends to behave better for difficult problems where the cost function is alternatively very flat and very steep. Argue for or against this conjecture.
(c) In support of your answer in (b), write a computer program to test the method with $\beta=0$ and $\beta>0$ with one-dimensional cost functions of the form

$$
f(x)=\frac{1}{2} x^{2}(1+\gamma \cos (x))
$$

where $\gamma \in(0,1)$, and

$$
f(x)=\frac{1}{2} \sum_{i=1}^{m}\left|z_{i}-\tanh \left(x y_{i}\right)\right|^{2}
$$

where $z_{i}$ and $y_{i}$ are given scalars.

# 1.4 NEWTON'S METHOD AND VARIATIONS 

In the last two sections we emphasized a basic tradeoff in gradient methods: implementation simplicity versus fast convergence. We have already discussed steepest descent, one of the simplest but also one of the slowest methods. We now consider its opposite extreme, Newton's method, which is arguably the most complex and also the fastest of the gradient methods (under appropriate conditions).

Newton's method consists of the iteration

$$
x^{k+1}=x^{k}-\alpha^{k}\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1} \nabla f\left(x^{k}\right)
$$

assuming that the Newton direction

$$
d^{k}=-\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1} \nabla f\left(x^{k}\right)
$$

---

is defined and is a direction of descent [i.e., $d_{k}{ }^{\prime} \nabla f\left(x^{k}\right)<0$ ]. As explained in the preceding section, one may view this iteration as a scaled version of steepest descent where the "optimal" scaling matrix $\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1}$ is used. It is worth mentioning in this connection that Newton's method is "scalefree", in the sense that it cannot be affected by a change in coordinate system as is true for steepest descent (see Exercise 4.1).

When the Armijo rule is used with initial stepsize $s=1$, then no reduction of the stepsize will be necessary near a nonsingular minimum (positive definite Hessian), as shown in Prop. 1.3.2. Thus, near convergence the method takes the form

$$
x^{k+1}=x^{k}-\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1} \nabla f\left(x^{k}\right)
$$

which will be referred to as the pure form of Newton's method. On the other hand, far from such a local minimum, the Hessian matrix may be singular or the Newton direction of Eq. (4.2) may not be a direction of descent because the Hessian $\nabla^{2} f\left(x^{k}\right)$ is not positive definite. Thus the analysis of Newton's method has two principal aspects:
(a) Local convergence, dealing with the behavior of the pure form of the method near a nonsingular local minimum.
(b) Global convergence, addressing the modifications that are necessary to ensure that the method is valid and is likely to converge to a local minimum when started far from all local minima.

We consider these issues in this section and we also discuss some variations of Newton's method, which are aimed at reducing the overhead for computing the Newton direction.

# Local Convergence 

The local convergence result for gradient methods (Prop. 1.2.5) together with the superlinear convergence result for Newton-like methods (Prop. 1.3.2) suggest that the pure form of Newton's method converges superlinearly when started close enough to a nonsingular local minimum. Results of this type hold for a more general form of Newton's method, that can be used to solve the system of $n$ equations with $n$ unknowns

$$
g(x)=0
$$

where $g: \Re^{n} \mapsto \Re^{n}$ is a continuously differentiable function. This method has the form

$$
x^{k+1}=x^{k}-\left(\nabla g\left(x^{k}\right)^{\prime}\right)^{-1} g\left(x^{k}\right)
$$

and for the special case where $g(x)$ is equal to the gradient $\nabla f(x)$, it yields the pure form of Eq. (4.3). [A continuously differentiable function $g: \Re^{n} \mapsto$

---

$\Re^{n}$ need not be equal to the gradient of some function. In particular, $g(x)=\nabla f(x)$ for some $f: \Re^{n} \mapsto \Re$, if and only if the $n \times n$ matrix $\nabla g(x)$ is symmetric for all $x([\mathrm{OrR} 70]$, p. 95$)$. Thus, the equation version of Newton's method (4.5) is more broadly applicable than the optimization version of Eq. (4.3).]

There is a simple argument that shows the fast convergence of Newton's method (4.5). Suppose that the method generates a sequence $\left\{x^{k}\right\}$ that converges to a vector $x^{*}$ such that $g\left(x^{*}\right)=0$ and $\nabla g\left(x^{*}\right)$ is invertible. Let us use Taylor's theorem to write

$$
0=g\left(x^{*}\right)=g\left(x^{k}\right)+\nabla g\left(x^{k}\right)^{\prime}\left(x^{*}-x^{k}\right)+o\left(\left\|x^{k}-x^{*}\right\|\right)
$$

By multiplying this relation with $\left(\nabla g\left(x^{k}\right)^{\prime}\right)^{-1}$ we have

$$
x^{k}-x^{*}-\left(\nabla g\left(x^{k}\right)^{\prime}\right)^{-1} g\left(x^{k}\right)=o\left(\left\|x^{k}-x^{*}\right\|\right)
$$

so for the pure Newton iteration $x^{k+1}=x^{k}-\left(\nabla g\left(x^{k}\right)^{\prime}\right)^{-1} g\left(x^{k}\right)$ we obtain

$$
x^{k+1}-x^{*}=o\left(\left\|x^{k}-x^{*}\right\|\right)
$$

or, for $x^{k} \neq x^{*}$,

$$
\lim _{k \rightarrow \infty} \frac{\left\|x^{k+1}-x^{*}\right\|}{\left\|x^{k}-x^{*}\right\|}=\lim _{k \rightarrow \infty} \frac{o\left(\left\|x^{k}-x^{*}\right\|\right)}{\left\|x^{k}-x^{*}\right\|}=0
$$

implying superlinear convergence. This argument can also be used to show convergence to $x^{*}$ if the initial vector $x^{0}$ is sufficiently close to $x^{*}$. The following proposition proves a more detailed result.

Proposition 1.4.1: Consider a function $g: \Re^{n} \mapsto \Re^{n}$, and a vector $x^{*}$ such that $g\left(x^{*}\right)=0$. For $\delta>0$, let $S_{\delta}$ denote the sphere $\{x \mid$ $\left.\left\|x-x^{*}\right\| \leq \delta\right\}$. Assume that $g$ is continuously differentiable within some sphere $S_{\bar{\delta}}$ and that $\nabla g\left(x^{*}\right)$ is invertible.
(a) There exists $\delta>0$ such that if $x^{0} \in S_{\delta}$, the sequence $\left\{x^{k}\right\}$ generated by the iteration

$$
x^{k+1}=x^{k}-\left(\nabla g\left(x^{k}\right)^{\prime}\right)^{-1} g\left(x^{k}\right)
$$

is defined, belongs to $S_{\delta}$, and converges to $x^{*}$. Furthermore, $\left\{\left\|x^{k}-x^{*}\right\|\right\}$ converges superlinearly.

---

(b) If for some $L>0, M>0, \delta>0$, and for all $x$ and $y$ in $S_{\delta}$,

$$
\|\nabla g(x)-\nabla g(y)\| \leq L\|x-y\|, \quad\left\|\left(\nabla g(x)^{\prime}\right)^{-1}\right\| \leq M
$$

then, if $x^{0} \in S_{\delta}$, we have

$$
\left\|x^{k+1}-x^{*}\right\| \leq \frac{L M}{2}\left\|x^{k}-x^{*}\right\|^{2}, \quad \forall k=0,1, \ldots
$$

so $\left\{\left\|x^{k}-x^{*}\right\|\right\}$ converges superlinearly with order at least two.

Proof: (a) Let $\delta>0$ be such that $\left(\nabla g(x)^{\prime}\right)^{-1}$ exists within $S_{\delta}$ and let $M>0$ be such that

$$
\left\|\left(\nabla g(x)^{\prime}\right)^{-1}\right\| \leq M, \quad \forall x \in S_{\delta}
$$

Assuming $x \in S_{\delta}$, and using the relation

$$
g\left(x^{k}\right)=\int_{0}^{1} \nabla g\left(x^{*}+t\left(x^{k}-x^{*}\right)\right)^{\prime} d t\left(x^{k}-x^{*}\right)
$$

we estimate $\left\|x^{k+1}-x^{*}\right\|$ as

$$
\begin{aligned}
& \left\|x^{k+1}-x^{*}\right\|=\left\|x^{k}-x^{*}-\left(\nabla g\left(x^{k}\right)^{\prime}\right)^{-1} g\left(x^{k}\right)\right\| \\
& \quad=\left\|\left(\nabla g\left(x^{k}\right)^{\prime}\right)^{-1}\left(\nabla g\left(x^{k}\right)^{\prime}\left(x^{k}-x^{*}\right)-g\left(x^{k}\right)\right)\right\| \\
& \quad=\left\|\left(\nabla g\left(x^{k}\right)^{\prime}\right)^{-1}\left(\nabla g\left(x^{k}\right)^{\prime}-\int_{0}^{1} \nabla g\left(x^{*}+t\left(x^{k}-x^{*}\right)\right)^{\prime} d t\right)\left(x^{k}-x^{*}\right)\right\| \\
& \quad=\left\|\left(\nabla g\left(x^{k}\right)^{\prime}\right)^{-1}\left(\int_{0}^{1}\left[\nabla g\left(x^{k}\right)^{\prime}-\nabla g\left(x^{*}+t\left(x^{k}-x^{*}\right)\right)^{\prime}\right] d t\right)\left(x^{k}-x^{*}\right)\right\| \\
& \quad \leq M\left(\int_{0}^{1}\left\|\nabla g\left(x^{k}\right)-\nabla g\left(x^{*}+t\left(x^{k}-x^{*}\right)\right)\right\| d t\right)\left\|x^{k}-x^{*}\right\|
\end{aligned}
$$

By continuity of $\nabla g$, we can take $\delta$ sufficiently small to ensure that the term under the integral sign is arbitrarily small. The convergence $x^{k} \rightarrow x^{*}$ and the superlinear convergence of $\left\|x^{k}-x^{*}\right\|$ follow.
(b) If the condition (4.6) holds, Eq. (4.7) yields

$$
\left\|x^{k+1}-x^{*}\right\| \leq M\left(\int_{0}^{1} L t\left\|x^{k}-x^{*}\right\| d t\right)\left\|x^{k}-x^{*}\right\|=\frac{L M}{2}\left\|x^{k}-x^{*}\right\|^{2}
$$

Q.E.D.

---

A related result is the following. Its proof is left for the reader.

Proposition 1.4.2: Under the assumptions of Prop. 1.4.1(a), given any $r>0$, there exists a $\delta>0$ such that if $\left\|x^{k}-x^{*}\right\|<\delta$, then

$$
\left\|x^{k+1}-x^{*}\right\| \leq r\left\|x^{k}-x^{*}\right\|, \quad\left\|g\left(x^{k+1}\right)\right\| \leq r\left\|g\left(x^{k}\right)\right\|
$$

Thus, once it gets "near" a solution $x^{*}$ where $\nabla g\left(x^{*}\right)$ is invertible, the pure form of Newton's method converges extremely fast, typically taking a handful of iterations to achieve very high solution accuracy; see Fig. 1.4.1. Unfortunately, it is typically difficult to predict whether a given starting point is sufficiently near to a solution for the fast convergence rate of Newton's method to become effective right away. Thus, in practice one can only expect that eventually the fast convergence rate of Newton's method will become effective. Figure 1.4.2 illustrates how the method can fail to converge when started far from a solution.

# Global Convergence 

Newton's method in its pure form for unconstrained minimization of $f$ has several serious drawbacks.
(a) The inverse $\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1}$ may fail to exist, in which case the method breaks down. This will happen, for example, in regions where $f$ is linear $\left(\nabla^{2} f=0\right)$.
(b) The pure form is not a descent method; it may happen that $f\left(x^{k+1}\right)>$ $f\left(x^{k}\right)$.
(c) The pure form tends to be attracted by local maxima just as much as it is attracted by local minima. It just tries to solve the system of equations $\nabla f(x)=0$.

For these reasons, it is necessary to modify the pure form of Newton's method to turn it into a reliable minimization algorithm. There are several schemes that accomplish this by converting the pure form into a gradient method with a gradient related direction sequence. Simultaneously the modifications are such that, near a nonsingular local minimum, the algorithm assumes the pure form of Newton's method (4.3) and achieves the attendant fast convergence rate.

A simple possibility is to replace the Newton direction by the steepest descent direction (possibly after diagonal scaling), whenever the Newton

---

![[chapter_21_p84_img22.jpeg]]

Figure 1.4.1. Fast convergence of Newton's method for solving the equation $e^{x}-1=0$.
![[chapter_21_p84_img23.jpeg]]

Figure 1.4.2. Divergence of Newton's method for solving an equation $g(x)=0$ of a single variable $x$, when the starting point is far from the solution. This phenomenon typically happens if $\|\nabla g(x)\|$ tends to decrease as $\|x\| \rightarrow \infty$.
direction is either not defined or is not a descent direction. $\dagger$ With proper
$\dagger$ Interestingly, this motivated the development of steepest descent by M. Augustin Cauchy. In his original paper [Cau47], Cauchy states as motivation for the steepest descent method its capability to obtain a close approximation to the solution, in which case " . . one can obtain new approximations very rapidly with

---

safeguards, such a method has appropriate convergence and asymptotic rate of convergence properties (see Exercise 4.3 and for a related method, see Exercise 4.4). However, its performance at the early iterations may be quite slow, whether the Newton direction or the steepest descent direction is used in these iterations.

Generally, no modification of Newton's method can be guaranteed to converge fast in the early iterations, but there are schemes that can use second derivative information effectively, even when the Hessian is not positive definite. These schemes are based on making diagonal modifications to the Hessian; that is, they obtain the direction $d^{k}$ by solving a system of the form

$$
\left(\nabla^{2} f\left(x^{k}\right)+\Delta^{k}\right) d^{k}=-\nabla f\left(x^{k}\right)
$$

whenever the Newton direction does not exist or is not a descent direction. Here $\Delta^{k}$ is a diagonal matrix such that

$$
\nabla^{2} f\left(x^{k}\right)+\Delta^{k}: \text { positive definite. }
$$

We outline some possibilities.

# Modified Cholesky Factorization* 

It can be shown that every positive definite matrix $Q$ has a unique factorization of the form

$$
Q=L L^{\prime}
$$

where $L$ is a lower triangular matrix; this is known as the Cholesky factorization of $Q$ (see Appendix D). Systems of equations of the form $Q x=b$ can be solved by first solving for $y$ the triangular system $L y=b$, and then by solving for $x$ the triangular system $L^{\prime} x=y$. These triangular systems can be solved easily [in $O\left(n^{2}\right)$ operations as opposed to general systems, which require $O\left(n^{3}\right)$ operations; see Appendix D]. Since calculation of the Newton direction involves solution of the system

$$
\nabla^{2} f\left(x^{k}\right) d^{k}=-\nabla f\left(x^{k}\right)
$$

it is natural to compute $d^{k}$ by attempting to form the Cholesky factorization of $\nabla^{2} f\left(x^{k}\right)$. During this process, one can detect whether $\nabla^{2} f\left(x^{k}\right)$ is either nonpositive definite or nearly singular, in which case some of the diagonal elements of $\nabla^{2} f\left(x^{k}\right)$ are suitably increased to ensure that the resulting matrix is positive definite. This is done sequentially during the factorization process, so in the end we obtain

$$
L^{k} L^{k^{\prime}}=\nabla^{2} f\left(x^{k}\right)+\Delta^{k}
$$

where $L^{k}$ is lower triangular and nonsingular, and $\Delta^{k}$ is diagonal.
the aid of the linear or Newton's method ..." (Note the attribution to Newton by Cauchy.)

---

As an illustration, consider the 2-dimensional case (for the general case, see Appendix D). Let

$$
\nabla^{2} f\left(x^{k}\right)=\left(\begin{array}{ll}
h_{11} & h_{12} \\
h_{21} & h_{22}
\end{array}\right)
$$

and let the desired factorization be of the form

$$
L L^{\prime}=\left(\begin{array}{cc}
\alpha & 0 \\
\gamma & \beta
\end{array}\right) \cdot\left(\begin{array}{cc}
\alpha & \gamma \\
0 & \beta
\end{array}\right)
$$

We choose $\alpha, \beta$, and $\gamma$, so that $\nabla^{2} f\left(x^{k}\right)=L L^{\prime}$ if $\nabla^{2} f\left(x^{k}\right)$ is positive definite, and we appropriately modify $h_{11}$ and $h_{22}$ otherwise. This determines the first diagonal element $\alpha$ according to the relation

$$
\alpha=\left\{\begin{array}{ll}
\sqrt{h_{11}} & \text { if } h_{11}>0 \\
\sqrt{h_{11}+\delta_{1}} & \text { otherwise }
\end{array}\right.
$$

where $\delta_{1}$ is such that $h_{11}+\delta_{1}>0$. Given $\alpha$, we can calculate $\gamma$ by equating the corresponding elements of $\nabla^{2} f\left(x^{k}\right)$ and $L L^{\prime}$. We obtain $\gamma \alpha=h_{12}$ or

$$
\gamma=\frac{h_{12}}{\alpha}
$$

We can now calculate the second diagonal element $\beta$ by equating the corresponding elements of $\nabla^{2} f\left(x^{k}\right)$ and $L L^{\prime}$, after appropriately modifying $h_{22}$ if necessary,

$$
\beta=\left\{\begin{array}{ll}
\sqrt{h_{22}-\gamma^{2}} & \text { if } h_{22}>\gamma^{2} \\
\sqrt{h_{22}-\gamma^{2}+\delta_{2}} & \text { otherwise }
\end{array}\right.
$$

where $\delta_{2}$ is such that $h_{22}-\gamma^{2}+\delta_{2}>0$. The method for choosing the increments $\delta_{1}$ and $\delta_{2}$ is largely heuristic. One possibility is discussed in Appendix D, which also describes more sophisticated versions of the above procedure where a positive increment is added to the diagonal elements of the Hessian even when the corresponding diagonal elements of the factorization are positive but very close to zero.

Given the $L^{k} L^{k^{\prime}}$ factorization, the direction $d^{k}$ is obtained by solving the system

$$
L^{k} L^{k^{\prime}} d^{k}=-\nabla f\left(x^{k}\right)
$$

The next iterate is

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}
$$

where $\alpha^{k}$ is chosen according to the Armijo rule or one of the other stepsize rules we have discussed.

To guarantee convergence, the increments added to the diagonal elements of the Hessian can be chosen so that $\left\{d^{k}\right\}$ is gradient related (cf. Prop. 1.2.1). Also, these increments can be chosen to be zero near a nonsingular local minimum. In particular, with proper safeguards, near such a point, the method becomes identical to the pure form of Newton's method and achieves the corresponding superlinear convergence rate (see Appendix D).

---

# Trust Region Methods* 

As explained in Section 1.2, the pure Newton step is obtained by minimizing over $d$ the second order Taylor series approximation of $f$ around $x^{k}$, given by

$$
f^{k}(d)=f\left(x^{k}\right)+\nabla f\left(x^{k}\right)^{\prime} d+\frac{1}{2} d^{\prime} \nabla^{2} f\left(x^{k}\right) d
$$

We know that $f^{k}(d)$ is a good approximation of $f\left(x^{k}+d\right)$ when $d$ is in a small neighborhood of zero, but the difficulty is that with unconstrained minimization of $f^{k}(d)$ one may obtain a step that lies outside this neighborhood. It therefore makes sense to consider a restricted Newton step $d^{k}$ obtained by minimizing $f^{k}(d)$ over a suitably small neighborhood of zero, called the trust region:

$$
d^{k}=\arg \min _{|d| \leq \gamma^{k}} f^{k}(d)
$$

where $\gamma^{k}$ is some positive scalar. It can be shown that the restricted Newton step $d^{k}$ also solves a system of the form $\left(\nabla^{2} f\left(x^{k}\right)+\Lambda^{k} I\right) d=-\nabla f\left(x^{k}\right)$, where $I$ is the identity matrix and $\Lambda^{k}$ is a nonnegative scalar (a Lagrange multiplier in the terminology of Chapter 3), so the preceding method of determining $d^{k}$ fits the general framework of using a correction of the Hessian matrix by a positive semidefinite matrix. An approximate solution of the constrained minimization problem of Eq. (4.8) can be obtained quickly using the fact that it has only one constraint (see [MoS83]).

An important observation here is that even if $\nabla^{2} f\left(x^{k}\right)$ is not positive definite or, more generally, even if the pure Newton direction is not a descent direction, the restricted Newton step $d^{k}$ improves the cost, provided $\nabla f\left(x^{k}\right) \neq$ 0 and $\gamma^{k}$ is sufficiently small. To see this, note that we have for all $d$ with $\|d\| \leq \gamma^{k}$

$$
f\left(x^{k}+d\right)=f^{k}(d)+o\left(\left(\gamma^{k}\right)^{2}\right)
$$

so that

$$
\begin{aligned}
f\left(x^{k}+d^{k}\right) & =f^{k}\left(d^{k}\right)+o\left(\left(\gamma^{k}\right)^{2}\right) \\
& =f\left(x^{k}\right)+\min _{\|d\| \leq \gamma^{k}}\left\{\nabla f\left(x^{k}\right)^{\prime} d+\frac{1}{2} d^{\prime} \nabla^{2} f\left(x^{k}\right) d\right\}+o\left(\left(\gamma^{k}\right)^{2}\right)
\end{aligned}
$$

Therefore, denoting

$$
\tilde{d}^{k}=-\frac{\nabla f\left(x^{k}\right)}{\left\|\nabla f\left(x^{k}\right)\right\|} \gamma^{k}
$$

we have

$$
\begin{aligned}
f\left(x^{k}+d^{k}\right) \leq & f\left(x^{k}\right)+\nabla f\left(x^{k}\right)^{\prime} \tilde{d}^{k}+\frac{1}{2} \tilde{d^{k}}^{\prime} \nabla^{2} f\left(x^{k}\right) \tilde{d}^{k}+o\left(\left(\gamma^{k}\right)^{2}\right) \\
= & f\left(x^{k}\right)-\gamma^{k}\left\|\nabla f\left(x^{k}\right)\right\|+\frac{\left(\gamma^{k}\right)^{2}}{2\left\|\nabla f\left(x^{k}\right)\right\|^{2}} \nabla f\left(x^{k}\right)^{\prime} \nabla^{2} f\left(x^{k}\right) \nabla f\left(x^{k}\right) \\
& +o\left(\left(\gamma^{k}\right)^{2}\right)
\end{aligned}
$$

For $\gamma^{k}$ sufficiently small, the negative term $-\gamma^{k}\left\|\nabla f\left(x^{k}\right)\right\|$ dominates the last two terms on the right-hand side above, showing that

$$
f\left(x^{k+1}\right)<f\left(x^{k}\right)
$$

---

It can be seen in fact from the preceding relations that a cost improvement is possible even when $\nabla f\left(x^{k}\right)=0$, provided $\gamma^{k}$ is sufficiently small and $f$ has a direction of negative curvature at $x^{k}$, that is, $\nabla^{2} f\left(x^{k}\right)$ is not positive semidefinite. Thus the preceding procedure will fail to improve the cost only if $\nabla f\left(x^{k}\right)=0$ and $\nabla^{2} f\left(x^{k}\right)$ is positive semidefinite, that is, $x^{k}$ satisfies the first and second order necessary conditions. In particular, one can typically make progress even if $x^{k}$ is a stationary point that is not a local minimum.

We are thus motivated to consider a method of the form

$$
x^{k+1}=x^{k}+d^{k}
$$

where $d^{k}$ is the restricted Newton step corresponding to a suitably chosen scalar $\gamma^{k}$ as per Eq. (4.8). Here, for a given $x^{k}, \gamma^{k}$ should be small enough so that there is cost improvement; one possibility is to start from an initial trial $\gamma^{k}$ and successively reduce $\gamma^{k}$ by a certain factor as many times as necessary until a cost reduction occurs $\left[f\left(x^{k+1}\right)<f\left(x^{k}\right)\right]$. The choice of the initial trial value for $\gamma^{k}$ is crucial here; if it is chosen too large, a large number of reductions may be necessary before a cost improvement occurs; if it is chosen too small the convergence rate may be poor. In particular, to maintain the superlinear convergence rate of Newton's method, as $x^{k}$ approaches a nonsingular local minimum, one should select the initial trial value of $\gamma^{k}$ sufficiently large so that the restricted Newton step and the pure Newton step coincide.

A reasonable way to adjust the initial trial value for $\gamma^{k}$ is to increase this value when the method appears to be progressing well and to decrease this value otherwise. One can measure progress by using the ratio of actual over predicted cost improvement [based on the approximation $f^{k}(d)$ ]

$$
r^{k}=\frac{f\left(x^{k}\right)-f\left(x^{k+1}\right)}{f\left(x^{k}\right)-f^{k}\left(d^{k}\right)}
$$

In particular, it makes sense to increase the initial trial value for $\gamma\left(\gamma^{k+1}>\gamma^{k}\right)$ if this ratio is close to or above unity, and decrease $\gamma$ otherwise. The following algorithm is a typical example of such a method. Given $x^{k}$ and an initial trial value $\gamma^{k}$, it determines $x^{k+1}$ and an initial trial value $\gamma^{k+1}$ by using two threshold values $\sigma_{1}, \sigma_{2}$ with $0<\sigma_{1} \leq \sigma_{2} \leq 1$ and two factors $\beta_{1}, \beta_{2}$ with $0<\beta_{1}<1<\beta_{2}$ (typical values are $\sigma_{1}=0.2, \sigma_{2}=0.8, \beta_{1}=0.25, \beta_{2}=2$ ).
Step 1: Find

$$
d^{k}=\arg \min _{||d|| \leq \gamma^{k}} f^{k}(d)
$$

If $f^{k}\left(d^{k}\right)=f\left(x^{k}\right)$ stop ( $x^{k}$ satisfies the first and second order necessary conditions for a local minimum); else go to Step 2.
Step 2: If $f\left(x^{k}+d^{k}\right)<f\left(x^{k}\right)$ set

$$
x^{k+1}=x^{k}+d^{k}
$$

calculate

$$
r^{k}=\frac{f\left(x^{k}\right)-f\left(x^{k+1}\right)}{f\left(x^{k}\right)-f^{k}\left(d^{k}\right)}
$$

---

and go to Step 3; else set $\gamma^{k}:=\beta_{1}\left\|d^{k}\right\|$ and go to Step 1.
Step 3: Set

$$
\gamma^{k+1}= \begin{cases}\beta_{1}\left\|d^{k}\right\| & \text { if } r^{k}<\sigma_{1} \\ \beta_{2} \gamma^{k} & \text { if } \sigma_{2} \leq r_{k} \text { and }\left\|d^{k}\right\|=\gamma^{k} \\ \gamma^{k} & \text { otherwise }\end{cases}
$$

Go to the next iteration.
Assuming that $f$ is twice continuously differentiable, it is possible to show that the above algorithm is convergent in the sense that if $\left\{x^{k}\right\}$ is a bounded sequence, there exists a limit point of $\left\{x^{k}\right\}$ that satisfies the first and the second order necessary conditions for optimality. Furthermore, if $\left\{x^{k}\right\}$ converges to a nonsingular local minimum $x^{*}$, then asymptotically, the method is identical to the pure form of Newton's method, thereby attaining a superlinear convergence rate; see the references given at the end of the chapter for proofs of these and other related results for trust region methods.

# Newton's Method with Periodic Reevaluation of the Hessian 

A variation of Newton's method is obtained if the Hessian matrix $\nabla^{2} f$ is recomputed every $p>1$ iterations rather than at every iteration. In particular, this method, in unmodified form, is given by

$$
x^{k+1}=x^{k}-\alpha^{k} D^{k} \nabla f\left(x^{k}\right)
$$

where

$$
D^{i p+j}=\left(\nabla^{2} f\left(x^{i p}\right)\right)^{-1}, \quad j=0,1, \ldots, p-1, i=0,1, \ldots
$$

The idea here is to save the computation and the inversion (or factorization) of the Hessian for the iterations where $j \neq 0$. This reduction in overhead is achieved at the expense of what is usually a small degradation in speed of convergence.

## Truncated Newton Methods*

We have so far implicitly assumed that the system $\nabla^{2} f\left(x^{k}\right) d^{k}=$ $-\nabla f\left(x^{k}\right)$ will be solved for the direction $d^{k}$ by Cholesky factorization or Gaussian elimination, which require a finite number of arithmetic operations $\left[O\left(n^{3}\right)\right]$. When the dimension $n$ is large, the calculation required for exact solution of the system may be prohibitive, and one may have to be satisfied with only an approximate solution. Such an approximation may be obtained by using an iterative method. This approach is often used for solving very large linear systems of equations, arising in the solution of partial differential equations, where an adequate approximation to the

---

solution can often be obtained by iterative methods quite fast, while the computation to find the exact solution can be overwhelming.

Generally, the solution in $d$ of any system of the form $H^{k} d=-\nabla f\left(x^{k}\right)$, where $H^{k}$ is a positive definite symmetric $n \times n$ matrix, is obtained by solving the quadratic optimization problem

$$
\begin{array}{ll}
\operatorname{minimize} & \frac{1}{2} d^{\prime} H^{k} d+\nabla f\left(x^{k}\right)^{\prime} d \\
\text { subject to } & d \in \Re^{n}
\end{array}
$$

whose cost function gradient is zero at $d$ if and only if $H^{k} d=-\nabla f\left(x^{k}\right)$. Suppose that an iterative descent method is used for solution and the starting point is $d^{0}=0$. Since the quadratic cost is reduced at each iteration and its value at the starting point is zero, we obtain after each iteration a vector $d^{k}$ satisfying $\frac{1}{2} d^{k^{\prime}} H^{k} d^{k}+\nabla f\left(x^{k}\right)^{\prime} d^{k}<0$, from which, using the positive definiteness of $H^{k}$,

$$
\nabla f\left(x^{k}\right)^{\prime} d^{k}<0
$$

Thus the approximate solution $d^{k}$ of the system $H^{k} d=-\nabla f\left(x^{k}\right)$, obtained after any positive number of iterations, is a descent direction. Possible iterative descent methods include the conjugate gradient method to be presented in Section 1.6 and the coordinate descent methods to be discussed in Section 1.8.

Conditions on the accuracy of the approximate solution $d^{k}$ that ensure linear or superlinear rate of convergence are given in Exercise 4.5. Generally, the superlinear convergence rate property of the method to a nonsingular local minimum is maintained if the approximate Newton directions $d^{k}$ satisfy

$$
\lim _{k \rightarrow \infty} \frac{\left\|\nabla^{2} f\left(x^{k}\right) d^{k}+\nabla f\left(x^{k}\right)\right\|}{\left\|\nabla f\left(x^{k}\right)\right\|}=0
$$

(compare with Prop. 1.3.2). Thus, for superlinear convergence rate, the norm of the residual error in solving the Newton system must become negligible relative to the gradient norm in the limit.

# E XERCISES 

4.1

Consider a linear invertible transformation of variables $x=S y$. Write Newton's method in the space of the variables $y$ and show that it generates the sequence $y^{k}=S^{-1} x^{k}$, where $\left\{x^{k}\right\}$ is the sequence generated by Newton's method in the space of the variables $x$.

---

4.2

Show Prop. 1.4.2. Hint: For the second relation, let $M(x)=\int_{0}^{1} \nabla g\left(x^{*}+\right.$ $\left.t\left(x-x^{*}\right)\right)^{\prime} d t$, so that $g(x)=M(x)\left(x-x^{*}\right)$. Argue that for some $\delta>0$ the eigenvalues of $M(x)^{\prime} M(x)$ lie between some positive scalars $\gamma$ and $I^{\prime}$ for all $x$ with $\left\|x-x^{*}\right\| \leq \delta$. Show that

$$
\gamma\left\|x-x^{*}\right\|^{2} \leq\|g(x)\|^{2} \leq \Gamma\left\|x-x^{*}\right\|^{2}, \quad \forall x \text { with }\left\|x-x^{*}\right\| \leq \delta
$$

# 4.3 (1st Combination with Steepest Descent) 

Consider the iteration $x^{k+1}=x^{k}+\alpha^{k} d^{k}$ where $\alpha^{k}$ is chosen by the Armijo rule with initial stepsize $s=1, \sigma \in(0,1 / 2)$, and $d^{k}$ is equal to

$$
d_{N}^{k}=-\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1} \nabla f\left(x^{k}\right)
$$

if $\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1}$ exists and

$$
\begin{gathered}
\nabla f\left(x^{k}\right)^{\prime} d_{N}^{k} \geq c_{1}\left\|\nabla f\left(x^{k}\right)\right\|^{p_{1}} \\
c_{2}\left\|\nabla f\left(x^{k}\right)\right\| \geq\left\|d_{N}^{k}\right\|^{p_{2}}
\end{gathered}
$$

while otherwise

$$
d^{k}=-D \nabla f\left(x^{k}\right)
$$

The matrix $D$ is positive definite and the scalars $c_{1}, c_{2}, p_{1}$, and $p_{2}$ satisfy

$$
c_{1}>0, \quad c_{2}>0, \quad p_{1}>2, \quad p_{2}>1
$$

Show that the sequence $\left\{d^{k}\right\}$ is gradient related. Furthermore if the method is started close enough to a nonsingular local minimum $x^{*}$, then $\left\{\left\|x^{k}-x^{*}\right\|\right\}$ converges to zero superlinearly.

## 4.4 (2nd Combination with Steepest Descent)

Assume that $f: \Re^{n} \mapsto \Re$ is twice continuously differentiablh with everywhere nonsingular (but not necessarily positive definite) Hessian. At a point $x^{k}$, let $d_{N}^{k}$ be the Newton direction and let $d_{S}^{k}=-D \nabla f\left(x^{k}\right)$ be a scaled steepest descent direction, where $D$ is a fixed positive definite symmetric matrix. Consider the method

$$
x^{k+1}=x^{k}+\alpha^{k}\left(\left(1-\alpha^{k}\right) d_{S}^{k}+\alpha^{k} d_{N}^{k}\right)
$$

where $\alpha^{k}=\beta^{m^{k}}$ and $m^{k}$ is the first nonnegative integer $m$ such that $f\left(x^{k}\right)-f\left(x^{k}+\beta^{m}\left(\left(1-\beta^{m}\right) d_{S}^{k}+\beta^{m} d_{N}\right)\right) \geq-\sigma \beta^{m} \nabla f\left(x^{k}\right)^{\prime}\left(\left(1-\beta^{m}\right) d_{S}^{k}+\beta^{m} d_{N}\right)$ and the scalars $\beta$ and $\sigma$ satisfy $0<\beta<1$ and $0<\sigma<1 / 2$. Show that the method is well defined in the sense that the stepsize $\alpha^{k}$ will be obtained after a finite number of trials. Furthermore, every limit point of $\left\{x^{k}\right\}$ is stationary and if $\left\{x^{k}\right\}$ converges to a nonsingular local minimum $x^{*}$, the rate of convergence of $\left\{\left\|x^{k}-x^{*}\right\|\right\}$ is superlinear.

---

4.5

Consider a truncated Newton method and assume that $\left\{x^{k}\right\}$ converges to a nonsingular local minimum $x^{*}$. Assume that the matrices $H^{k}$ and the directions $d^{k}$ satisfy

$$
\lim _{k \rightarrow \infty}\left\|H^{k}-\nabla^{2} f\left(x^{k}\right)\right\|=0, \quad \lim _{k \rightarrow \infty} \frac{\left\|H^{k} d^{k}+\nabla f\left(x^{k}\right)\right\|}{\left\|\nabla f\left(x^{k}\right)\right\|}=0
$$

Show that $\left\{\left\|x^{k}-x^{*}\right\|\right\}$ converges superlinearly.

# 4.6 

Apply Newton's method to minimization of the function $f(x)=\|x\|^{3}$ and show that it converges linearly to $x^{*}=0$. Explain this fact in light of Prop. 1.4.1.
4.7

Apply Newton's method with the trust region implementation to a positive definite quadratic function. Show that the method terminates in a finite number of iterations.

### 1.5 LEAST SQUARES PROBLEMS

In this section we consider methods for solving least squares problems of the form

$$
\begin{aligned}
& \operatorname{minimize} \quad f(x)=\frac{1}{2}\|g(x)\|^{2}=\frac{1}{2} \sum_{i=1}^{m}\left\|g_{i}(x)\right\|^{2} \\
& \text { subject to } \quad x \in \Re^{n}
\end{aligned}
$$

where $g$ is a continuously differentiable function with component functions $g_{1}, \ldots, g_{m}$, where $g_{i}: \Re^{n} \rightarrow \Re^{r_{i}}$. Usually $r_{i}=1$, but it is sometimes notationally convenient to consider the more general case.

Least squares problems are very common in practice. A principal case arises when $g$ consists of $n$ scalar-valued functions and we want to solve the system of $n$ equations with $n$ unknowns $g(x)=0$. We can formulate this as the least squares optimization problem (5.1) [ $x^{*}$ solves the system $g(x)=0$ if and only if it minimizes $\frac{1}{2}\|g(x)\|^{2}$ and the optimal value is zero]. Here are some other examples:

---

# Example 5.1 (Model Construction - Curve Fitting) 

Suppose that we want to estimate $n$ parameters of a mathematical model so that it fits well a physical system, based on a set of measurements. In particular, we hypothesize an approximate relation of the form

$$
z=h(x, y)
$$

where $h$ is a known function representing the model and
$x \in \Re^{n}$ is a vector of unknown parameters,
$z \in \Re^{r}$ is the model's output,
$y \in \Re^{p}$ is the model's input.
Given a set of $m$ input-output data pairs $\left(y_{1}, z_{1}\right), \ldots,\left(y_{m}, z_{m}\right)$ from measurements of the physical system that we try to model, we want to find the vector of parameters $x$ that matches best the data in the sense that it minimizes the sum of squared errors

$$
\frac{1}{2} \sum_{i=1}^{m}\left\|z_{i}-h\left(x, y_{i}\right)\right\|^{2}
$$

For example, to fit the data pairs by a cubic polynomial approximation, we would choose

$$
h(x, y)=x_{3} y^{3}+x_{2} y^{2}+x_{1} y+x_{0}
$$

where $x=\left(x_{0}, x_{1}, x_{2}, x_{3}\right)$ is the vector of unknown coefficients of the cubic polynomial.

The next two examples are really special cases of the preceding one.

## Example 5.2 (Dynamic System Identification)

A common model for a single input-single output dynamic system is to relate the input sequence $\left\{y_{k}\right\}$ to the output sequence $\left\{z_{k}\right\}$ by a linear equation of the form

$$
\sum_{j=0}^{n} \alpha_{j} z_{k-j}=\sum_{j=0}^{n} \beta_{j} y_{k-j}
$$

Given a record of inputs and outputs $y_{1}, z_{1}, \ldots, y_{m}, z_{m}$ from the true system, we would like to find a set of parameters $\left\{\alpha_{j}, \beta_{j} \mid j=0,1, \ldots, n\right\}$ that matches this record best in the sense that it minimizes

$$
\sum_{k=n}^{m}\left(\sum_{j=0}^{n} \alpha_{j} z_{k-j}-\sum_{j=0}^{n} \beta_{j} y_{k-j}\right)^{2}
$$

This is a least-squares problem.

---

# Example 5.3 (Neural Networks) 

A least squares modeling problem that has received a lot of attention is provided by neural networks. Here the model is specified by a multistage system, also called a multilayer perceptron. The $k$ th stage consists of $n_{k}$ activation units, each of which is a single input-single output mapping of a given form $\phi: \Re \mapsto \Re$ to be described shortly. The output of the $j$ th activation unit of the $(k+1)$ st stage is denoted by $x_{k+1}^{j}$ and the input is a linear function of the output vector $x_{k}=\left(x_{k}^{1}, \ldots, x_{k}^{n_{k}}\right)$ of the $k$ th stage. Thus

$$
x_{k+1}^{j}=\phi\left(u_{k}^{0 j}+\sum_{s=1}^{n_{k}} x_{k}^{s} u_{k}^{s j}\right), \quad j=1, \ldots, n_{k+1}
$$

where the coefficients $u_{k}^{s j}$ (also called weights) are to be determined.
Suppose that the multilayer perceptron has $N$ stages, and let $u$ denote the vector of the weights of all the stages:

$$
u=\left\{u_{k}^{s j} \mid k=0, \ldots, N-1, s=0, \ldots, n_{k}, j=1, \ldots, n_{k+1}\right\}
$$

Then, for a given vector $u$ of weights, an input vector $x_{0}$ to the first stage produces a unique output vector $x_{N}$ from the $N$ th stage via Eq. (5.2). Thus, we may view the multilayer perceptron as a mapping $h$ that is parameterized by $u$ and transforms the input vector $x_{0}$ into an output vector of the form $x_{N}=h\left(u, x_{0}\right)$. Suppose that we have $m$ sample input-output pairs $\left(y_{1}, z_{1}\right), \ldots,\left(y_{m}, z_{m}\right)$ from a physical system that we are trying to model. Then, by selecting $u$ appropriately, we can try to match the mapping of the multilayer perceptron with the mapping of the physical system. A common way to do this is to minimize over $u$ the sum of squared errors

$$
\frac{1}{2} \sum_{i=1}^{m}\left\|z_{i}-h\left(u, y_{i}\right)\right\|^{2}
$$

In the terminology of neural network theory, the process of finding the optimal weights is known as training the network.

Common examples of activation units are functions such as

$$
\begin{gathered}
\phi(\xi)=\frac{1}{1+e^{-\xi}}, \quad \text { (sigmoidal function) } \\
\phi(\xi)=\frac{e^{\xi}-e^{-\xi}}{e^{\xi}+e^{-\xi}}, \quad \text { (hyperbolic tangent function) }
\end{gathered}
$$

whose gradients are zero as the argument $\xi$ approaches $-\infty$ and $\infty$. For these functions $\phi$, it is possible to show that with a sufficient number of activation units and a number of stages $N \geq 2$, a multilayer perceptron can approximate arbitrarily closely very complex input-output maps; see [Cyb89].

Neural network training problems can be quite challenging. Their cost function is typically nonconvex and involves multiple local minima. For large

---

![[chapter_21_p95_img24.jpeg]]

Figure 1.5.1. Three-dimensional plot of a least squares cost function

$$
\frac{1}{2} \sum_{i=1}^{5}\left(z_{i}-\phi\left(u_{1} y_{i}+u_{0}\right)\right)^{2}
$$

for a neural network training problem where there are only two weights $u_{0}$ and $u_{1}$, five data pairs, and $\phi$ is the hyperbolic tangent function. The data of the problem are given in Exercise 5.3. The cost function tends to a constant as $u$ is changed along rays of the form $r \bar{u}$, where $r>0$ and $\bar{u}$ is a fixed vector.
values of the weights $u_{k}^{17}$, the cost becomes "flat". In fact, as illustrated in Fig. 1.5.1, the cost function tends to a constant as $u$ is changed along rays of the form $r \bar{u}$, where $r>0$ and $\bar{u}$ is a fixed vector. For $u$ near the origin, the cost function can be quite complicated alternately involving flat and steep regions.

The next example deals with an important context where neural networks are often used:

# Example 5.4 (Pattern Classification) 

Consider the problem of classifying objects based on the values of their char-

---

acteristics. (We use the term "object" generically; in some contexts, the classification may relate to persons or situations.) Each object is presented to us with a vector $y$ of features, and we wish to classify it in one of $s$ categories $1, \ldots, s$. For example, the vector $y$ may represent the results of a collection of tests on a medical patient, and we may wish to classify the patient as being healthy or as having one of several types of illnesses.

A classical pattern classification approach is to assume that we know the probabilities $p(j \mid y)$ of an object with feature vector $y$ being of category $j$, where $j=1, \ldots, s$. Then we may associate an object with feature vector $y$ with the category $j^{*}(y)$ having maximum posterior probability, that is,

$$
j^{*}(y)=\arg \max _{j=1, \ldots, s} p(j \mid y)
$$

Suppose now that the probabilities $p(j \mid y)$ are unknown, but instead we have a sample consisting of $m$ object-category pairs. Then we may try to estimate $p(j \mid y)$ based on the fact that, out of all functions $f_{j}(y)$ of $y, p(j \mid y)$ is the one that minimizes the expected value of $\left(z_{j}-f_{j}(y)\right)^{2}$, where

$$
z_{j}= \begin{cases}1 & \text { if } y \text { is of category } j \\ 0 & \text { otherwise }\end{cases}
$$

In particular, let $y_{i}$ denote the feature vector of the $i$ th object. For each category $j=1, \ldots, s$, we estimate the probability $p(j \mid y)$ as a function of $y$, by a function $h_{j}\left(x_{j}, y\right)$ that is parameterized by a vector $x_{j}$. The function $h_{j}$ could be provided for example by a neural network (cf. Example 5.3). Then, we can obtain $x_{j}$ by minimizing the least squares function

$$
\frac{1}{2} \sum_{i=1}^{m}\left(z_{j}^{i}-h_{j}\left(x_{j}, y_{i}\right)\right)^{2}
$$

where

$$
z_{j}^{i}= \begin{cases}1 & \text { if } y_{i} \text { is of category } j \\ 0 & \text { otherwise }\end{cases}
$$

This minimization approximates the minimization of the expected value of $\left(z_{j}-f_{j}(y)\right)^{2}$. Once the optimal parameter vectors $x_{j}^{*}, j=1, \ldots, s$, have been obtained, we may use them to classify a new object with feature vector $y$ according to the rule

$$
\text { Estimated Object Category }=\arg \max _{j=1, \ldots, s} h_{j}\left(x_{j}^{*}, y\right)
$$

which approximates the maximum posterior probability rule (5.3).
For the simpler case where there are just two categories, say $A$ and $B$, a similar formulation is to hypothesize a relation of the following form between feature vector $y$ and category of an object:

$$
\text { Object Category }= \begin{cases}A & \text { if } h(x, y)=1 \\ B & \text { if } h(x, y)=-\mathbf{1}\end{cases}
$$

---

where $h$ is a given function and $x$ is an unknown vector of parameters. Given a set of $m$ data pairs $\left(z_{1}, y_{1}\right), \ldots,\left(z_{m}, y_{m}\right)$ of representative objects of known category, where $y_{i}$ is the feature vector of the $i$ th object, and

$$
z_{i}= \begin{cases}1 & \text { if } y \text { is of category } A \\ -1 & \text { if } y \text { is of category } B\end{cases}
$$

we obtain $x$ by minimizing the least squares function

$$
\frac{1}{2} \sum_{i=1}^{m}\left(z_{i}-h\left(x, y_{i}\right)\right)^{2}
$$

The optimal parameter vector $x^{*}$ is used to classify a new object with feature vector $y$ according to the rule

$$
\text { Estimated Object Category }= \begin{cases}A & \text { if } h\left(x^{*}, y\right)>0 \\ B & \text { if } h\left(x^{*}, y\right)<0\end{cases}
$$

There are several other variations on the above theme, for which we refer to the specialized literature.

# 1.5.1 The Gauss-Newton Method 

Let us consider now specialized methods for minimizing the least squares cost $(1 / 2) \| g(x) \|^{2}$, starting with the most commonly used method, the Gauss-Newton method. Given a point $x^{k}$, the pure form of the GaussNewton iteration is based on linearizing $g$ to obtain

$$
\tilde{g}\left(x, x^{k}\right)=g\left(x^{k}\right)+\nabla g\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)
$$

and then minimizing the norm of the linearized function $\tilde{g}$ :

$$
\begin{aligned}
x^{k+1}= & \arg \min _{x \in \Re^{n}} \frac{1}{2}\left\|\tilde{g}\left(x, x^{k}\right)\right\|^{2} \\
= & \arg \min _{x \in \Re^{n}} \frac{1}{2}\left\{\left\|g\left(x^{k}\right)\right\|^{2}+2\left(x-x^{k}\right)^{\prime} \nabla g\left(x^{k}\right) g\left(x^{k}\right)\right. \\
& \left.\quad+\left(x-x^{k}\right)^{\prime} \nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)\right\}
\end{aligned}
$$

Assuming that the $n \times n$ matrix $\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}$ is invertible, the above quadratic minimization yields

$$
x^{k+1}=x^{k}-\left(\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}\right)^{-1} \nabla g\left(x^{k}\right) g\left(x^{k}\right)
$$

Note that if $g$ is already a linear function, we have $\|g(x)\|^{2}=\left\|\tilde{g}\left(x, x^{k}\right)\right\|^{2}$, and the method converges in a single iteration. Note also that the direction

$$
-\left(\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}\right)^{-1} \nabla g\left(x^{k}\right) g\left(x^{k}\right)
$$

---

used in the above iteration is a descent direction since $\nabla g\left(x^{k}\right) g\left(x^{k}\right)$ is the gradient at $x^{k}$ of the least squares cost function $(1 / 2)\|g(x)\|^{2}$ and $\left(\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}\right)^{-1}$ is a positive definite matrix.

To ensure descent, and also to deal with the case where the matrix $\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}$ is singular (as well as enhance convergence when this matrix is nearly singular), the method is often implemented in the modified form

$$
x^{k+1}=x^{k}-\alpha^{k}\left(\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}+\Delta^{k}\right)^{-1} \nabla g\left(x^{k}\right) g\left(x^{k}\right)
$$

where $\alpha^{k}$ is a stepsize chosen by one of the stepsize rules that we have discussed, and $\Delta^{k}$ is a diagonal matrix such that

$$
\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}+\Delta^{k}: \text { positive definite. }
$$

For example, $\Delta^{k}$ may be chosen in accordance with the Cholesky factorization scheme outlined in Section 1.4. An early proposal, known as the Levenberg-Marquardt method, is to choose $\Delta^{k}$ to be a positive multiple of the identity matrix. With these choices of $\Delta^{k}$, it can be seen that the directions used by the method are gradient related, and the convergence results of Section 1.2 apply.

# Relation to Newton's Method 

The Gauss-Newton method bears a close relation to Newton's method. In particular, assuming each $g_{i}$ is a scalar function, the Hessian of the cost $(1 / 2)\|g(x)\|^{2}$ is

$$
\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}+\sum_{i=1}^{m} \nabla^{2} g_{i}\left(x^{k}\right) g_{i}\left(x^{k}\right)
$$

so it is seen that the Gauss-Newton iterations (5.4) and (5.5) are approximate versions of their Newton counterparts, where the second order term

$$
\sum_{i=1}^{m} \nabla^{2} g_{i}\left(x^{k}\right) g_{i}\left(x^{k}\right)
$$

is neglected. Thus, in the Gauss-Newton method, we save the computation of this term at the expense of some deterioration in the convergence rate. If, however, the neglected term (5.7) is relatively small near a solution, the convergence rate of the Gauss-Newton method is satisfactory. This is often true in many applications such as for example when $g$ is nearly linear, and also when the components $g_{i}(x)$ are small near the solution. In the case when $m=n$ and we try to solve the system $g(x)=0$, the neglected term is zero at a solution. In this case, assuming $\nabla g\left(x^{k}\right)$ is invertible, we have

$$
\left(\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}\right)^{-1} \nabla g\left(x^{k}\right) g\left(x^{k}\right)=\left(\nabla g\left(x^{k}\right)^{\prime}\right)^{-1} g\left(x^{k}\right)
$$

---

and the pure form of the Gauss-Newton method (5.4) takes the form

$$
x^{k+1}=x^{k}-\left(\nabla g\left(x^{k}\right)^{\prime}\right)^{-1} g\left(x^{k}\right)
$$

which is identical to Newton's method for solving the system $g(x)=0$ [rather than Newton's method for minimizing $\|g(x)\|^{2}$ ]. Thus, the convergence rate is typically superlinear in this case, as discussed in Section 1.4 .

# 1.5.2 Incremental Gradient Methods* 

Let us return to the model construction Example 5.1, where we want to find a vector $x \in \Re^{n}$ of model parameters based on data obtained from a physical system. Each component $g_{i}$ in the least squares formulation is referred to as a data block, and the entire function $g=\left(g_{1}, \ldots, g_{m}\right)$ is referred to as the data set.

In many problems of interest where there are many data blocks, the Gauss-Newton method may be ineffective, because the size of the data set makes each iteration very costly. For such problems it may be more attractive to use an incremental method that does not wait to process the entire data set before updating $x$ : instead, the method cycles through the data blocks in sequence and updates the estimate of $x$ after each data block is processed.

For example, given $x^{k}$ we may obtain $x^{k+1}$ at the last step of the following algorithm

$$
\psi_{i}=\psi_{i-1}-\alpha^{k} h_{i}, \quad i=1, \ldots, m
$$

where

$$
\psi_{0}=x^{k}
$$

$\alpha^{k}>0$ is a stepsize, and the direction $h_{i}$ is the gradient of the $i$ th data block,

$$
h_{i}=\nabla g_{i}\left(\psi_{i-1}\right) g_{i}\left(\psi_{i-1}\right)
$$

This method, assuming the same stepsize $\alpha^{k}$ is used for all $i$, can be written as

$$
x^{k+1}=x^{k}-\alpha^{k} \sum_{i=1}^{m} \nabla g_{i}\left(\psi_{i-1}\right) g_{i}\left(\psi_{i-1}\right)
$$

It may be viewed as an incremental version of the steepest descent method, which is

$$
x^{k+1}=x^{k}-\alpha^{k} \nabla f\left(x^{k}\right)=x^{k}-\alpha^{k} \sum_{i=1}^{m} \nabla g_{i}\left(x^{k}\right) g_{i}\left(x^{k}\right)
$$

---

The incremental algorithm (5.10) has been used extensively for training of neural networks (see Example 5.3). We refer to it as an incremental gradient method. Its advantage is that it may converge much faster than the corresponding steepest descent method, particularly when far from the eventual limit. This type of behavior is most vividly illustrated in the case where the data blocks are linear and the vector $x$ is one-dimensional.

# Example 5.5 

Assume that $x$ is a scalar, and that the least squares problem has the form

$$
\begin{aligned}
& \operatorname{minimize} \quad f(x)=\frac{1}{2} \sum_{i=1}^{m}\left(a_{i} x-b_{i}\right)^{2} \\
& \text { subject to } \quad x \in \Re
\end{aligned}
$$

where $a_{i}$ and $b_{i}$ are given scalars with $a_{i} \neq 0$ for all $i$. The minimum of each of the squared data blocks

$$
f_{i}(x)=\frac{1}{2}\left(a_{i} x-b_{i}\right)^{2}
$$

is

$$
x_{i}^{*}=\frac{b_{i}}{a_{i}}
$$

while the minimum of the least squares cost function $f$ is

$$
x^{*}=\frac{\sum_{i=1}^{m} a_{i} b_{i}}{\sum_{i=1}^{m} a_{i}^{2}}
$$

It can be seen that $x^{*}$ lies within the range of the data block minima

$$
R=\left[\min _{i} x_{i}^{*}, \max _{i} x_{i}^{*}\right]
$$

and that for all $x$ outside the range $R$, the gradient

$$
\nabla f_{i}(x)=a_{i}\left(a_{i} x-b_{i}\right)
$$

has the same sign as $\nabla f(x)$ (see Fig. 1.5.2). As a result, when outside the region $R$, the incremental gradient method approaches $x^{*}$ at each step [cf. Eq. (5.8)]

$$
\psi_{i}=\psi_{i-1}-\alpha^{k} a_{i}\left(a_{i} \psi_{i-1}-b_{i}\right)
$$

provided the stepsize $\alpha^{k}$ is small enough. In fact it is sufficient that

$$
\alpha^{k} \leq \min _{i} \frac{1}{a_{i}^{2}}
$$

However, for $x$ inside the region $R$, the $i$ th step of a cycle of the incremental gradient method need not make progress. It will approach $x^{*}$ (for small

---

enough stepsize $\alpha^{k}$ ) only if the current point $\psi_{i-1}$ does not lie in the interval connecting $x_{i}^{*}$ and $x^{*}$. This induces an oscillatory behavior within the region $R$, and as a result, the incremental gradient method will typically not converge to $x^{*}$ unless $\alpha^{k} \rightarrow 0$. By contrast, it can be shown that the steepest descent method, which takes the form

$$
x^{k+1}=x^{k}-\alpha^{k} \sum_{i=1}^{m} a_{i}\left(a_{i} x^{k}-b_{i}\right)
$$

converges to $x^{*}$ for any constant stepsize satisfying

$$
\alpha^{k} \leq \frac{1}{\sum_{i=1}^{m} a_{i}^{2}}
$$

However, unless the stepsize choice is particularly favorable, for $x$ outside the region $R$, a full iteration of steepest descent need not make more progress towards the solution than a single step of the incremental gradient method. In other words, with comparably intelligent stepsize choices, far from the solution (outside $R$ ), a single pass through the entire data set by incremental gradient is roughly as effective as $m$ passes through the data set by steepest descent.
![[chapter_21_p101_img25.jpeg]]

Figure 1.5.2. Illustrating the advantage of incrementalism when far from the optimal solution. The $i$ th step in an incremental gradient cycle is a gradient step for minimizing $\left(a_{i} x-b_{i}\right)^{2}$, so if $x$ lies outside the region of data block minima

$$
R=\left[\min _{i} x_{i}^{*}, \max _{i} x_{i}^{*}\right]
$$

and the stepsize is small enough, progress towards the solution $x^{*}$ is made.
The preceding example relies on $x$ being one-dimensional, but in many multidimensional problems the same qualitative behavior can be observed.

---

In particular, a pass through the $i$ th data block $g_{i}$ by the incremental gradient method can make progress towards the solution in the region where the data block gradient $\nabla g_{i}\left(\psi_{i-1}\right) g_{i}\left(\psi_{i-1}\right)$ makes an angle less than 90 degrees with the cost function gradient $\nabla f\left(\psi_{i-1}\right)$. If the data blocks $g_{i}$ are not "too dissimilar", this is likely to happen in a region of points that are not too close to the optimal solution set; see also Exercise 5.5.

The choice of the stepsize $\alpha^{k}$ plays an important role in the performance of incremental gradient methods. On close examination, it turns out that the direction used by the method differs from the gradient direction by an error that can be much larger than the gradient, and for this reason a diminishing stepsize is essential for convergence to a stationary point of $f$. However, it turns out that a peculiar form of convergence also typically occurs for a constant but sufficiently small stepsize. In this case, the iterates converge to a "limit cycle", whereby the $i$ th iterates $\psi_{i}$ within the cycles converge to a different limit that the $j$ th iterates $\psi_{j}$ for $i \neq j$. The sequence $\left\{x^{k}\right\}$ that consists of the iterates obtained at the end of cycles converges, except that the limit obtained need not be a stationary point of $f$ (see Exercise 5.2 and the following analysis). The limit tends to be close to a stationary point when the constant stepsize is small (see the following Prop. 1.5.1). In practice, it is common to use a constant stepsize for a (possibly prespecified) number of iterations, then decrease the stepsize by a certain factor, and repeat, up to the point where the stepsize reaches a prespecified minimum. An alternative possibility is to use a stepsize rule of the form

$$
\alpha^{k}=\min \left\{\gamma, \frac{\gamma_{1}}{k+\gamma_{2}}\right\}
$$

where $\gamma, \gamma_{1}$, and $\gamma_{2}$ are some positive scalars. There are also variants of the incremental gradient method that use a constant stepsize throughout, and generically converge to a stationary point of $f$ at a linear rate. The degree of incrementalism in these methods gradually diminishes as the method progresses (see Exercise 5.4 and [Ber95b]).

We note a popular modification of the incremental gradient method, which uses the update

$$
\psi_{i}=\psi_{i-1}-\alpha^{k} h_{i}+\beta\left(\psi_{i-1}-\psi_{i-2}\right), \quad i=1, \ldots, m
$$

where $\beta \in[0,1)$, in place of the update $\psi_{i}=\psi_{i-1}-\alpha^{k} h_{i}$. This method may be viewed as an incremental version of the heavy ball method given in Exercise 3.9 of Section 1.3. The term $\left(\psi_{i-1}-\psi_{i-2}\right)$ is known as a momentum term, and is often helpful in dealing with the peculiar features of the cost functions of neural network training problems.

Another popular technique for incremental methods is to reshuffle randomly the order of the data blocks after each cycle through the data set. A related alternative is to select randomly the data block from the data set at each iteration. If we take the view that an incremental method

---

is basically a gradient method with errors, we see that randomization of the order of the data blocks tends to randomize the size of the errors, and it appears that this tends to improve the convergence properties of the method. This is also supported by analysis of gradient methods with random errors (see e.g. [BeT89], [KuC78], [PoT73a], [Pol87], [TBA86]).

Generally, it may be said that while incremental methods are used widely in practice, particularly in neural network training problems, their effective use often requires skill, insight into the problem's structure, and trial and error.

# Convergence Analysis of Incremental Gradient Methods 

The incremental gradient method of Eq. (5.10) admits a similar analysis to gradient methods with errors. In Prop. 1.5.1 below we consider the case where the data blocks are linear. The main idea of the proof is that the incremental gradient method can be viewed as the regular steepest descent iteration where the gradient is perturbed by an error term that is proportional to the stepsize [see Eq. (5.26) below]. The qualitative behavior of the method then is as discussed in Section 1.2. The following lemma will be used in the proof.

Lemma 1.5.1: Suppose that $\left\{e^{k}\right\}$ and $\left\{\gamma^{k}\right\}$ are nonnegative sequences, and $c$ is a positive constant such that

$$
e^{k+1} \leq\left(1-\gamma^{k}\right) e^{k}+c\left(\gamma^{k}\right)^{2}, \quad \gamma^{k} \leq 1, \quad \forall k=0,1, \ldots
$$

and

$$
\gamma^{k} \rightarrow 0, \quad \sum_{k=0}^{\infty} \gamma^{k}=\infty
$$

Then $e^{k} \rightarrow 0$.

Proof: We first show that given any $\epsilon>0$, we have $e^{k}<\epsilon$ for infinitely many $k$. Indeed, assuming this were not so and letting $\bar{k}$ be such that $e^{k} \geq \epsilon$ and $c \gamma^{k} \leq \epsilon / 2$ for all $k \geq \bar{k}$, we would have for all $k \geq \bar{k}$

$$
e^{k+1} \leq e^{k}-\gamma^{k} e^{k}+c\left(\gamma^{k}\right)^{2} \leq e^{k}-\gamma^{k} \epsilon+\gamma^{k} \epsilon / 2=e^{k}-\gamma^{k} \epsilon / 2
$$

Therefore, for all $m \geq \bar{k}$ we would have

$$
e^{m+1} \leq e^{\bar{k}}-(\epsilon / 2) \sum_{k=\bar{k}}^{m} \gamma^{k}
$$

This contradicts the nonnegativity of $\left\{e^{k}\right\}$ and the assumption $\sum_{k=0}^{\infty} \gamma^{k}=$ $\infty$.

---

Thus, given any $\epsilon>0$, there exists $\bar{k}$ such that $c \gamma^{k}<\epsilon$ for all $k \geq \bar{k}$ and $\epsilon^{\bar{k}}<\epsilon$. We then have

$$
e^{\bar{k}+1} \leq\left(1-\gamma^{k}\right) e^{\bar{k}}+c\left(\gamma^{k}\right)^{2}<\left(1-\gamma^{k}\right) \epsilon+\gamma^{k} \epsilon=\epsilon
$$

By repeating this argument, we obtain $e^{k}<\epsilon$ for all $k \geq \bar{k}$. Since $\epsilon$ can be arbitrarily small, it follows that $e^{k} \rightarrow 0$. Q.E.D.

Proposition 1.5.1: Consider the case of linear data blocks,

$$
g_{i}(x)=z_{i}-C_{i} x, \quad i=1, \ldots, m
$$

and the incremental gradient method

$$
x^{k+1}=x^{k}+\alpha^{k} \sum_{i=1}^{m} C_{i}^{\prime}\left(z_{i}-C_{i} \psi_{i-1}\right)
$$

where $\psi_{0}=x^{k}$ and

$$
\psi_{i}=\psi_{i-1}+\alpha^{k} C_{i}^{\prime}\left(z_{i}-C_{i} \psi_{i-1}\right), \quad i=1, \ldots, m
$$

Assume that $\sum_{i=1}^{m} C_{i}^{\prime} C_{i}$ is a positive definite matrix and let $x^{*}$ be the optimal solution of the corresponding least squares problem. Then:
(a) There exists $\bar{\alpha}>0$ such that if $\alpha^{k}$ is equal to some constant $\alpha \in$ $(0, \bar{\alpha}]$ for all $k,\left\{x^{k}\right\}$ converges to some vector $x(\alpha)$. Furthermore, we have $\lim _{\alpha \rightarrow 0} x(\alpha)=x^{*}$.
(b) If $\alpha^{k}>0$ for all $k$, and

$$
\alpha^{k} \rightarrow 0, \quad \sum_{k=0}^{\infty} \alpha^{k}=\infty
$$

then $\left\{x^{k}\right\}$ converges to $x^{*}$.

Proof: (a) We first show by induction that for all $i=1, \ldots, m$, the vectors $\psi_{i}$ of Eq. (5.14) have the form

$$
\psi_{i}=x^{k}+\alpha \sum_{j=1}^{i} C_{j}^{\prime}\left(z_{j}-C_{j} x^{k}\right)+\sum_{j=1}^{i-1} \alpha^{j+1}\left(\Phi_{i j} x^{k}+\phi_{i j}\right)
$$

where $\Phi_{i j}$ and $\phi_{i j}$ are some matrices and vectors, respectively, which are independent of $k$. Indeed, this relation holds by definition for $i=1$ [the

---

last term in Eq. (5.15) is zero when $i=1$ ]. Suppose that it holds for $i=p$. We have, using the induction hypothesis,

$$
\begin{aligned}
\psi_{p+1}= & \psi_{p}+\alpha C_{p+1}^{\prime}\left(z_{p+1}-C_{p+1} \psi_{p}\right) \\
= & x^{k}+\alpha \sum_{j=1}^{p} C_{j}^{\prime}\left(z_{j}-C_{j} x^{k}\right)+\sum_{j=1}^{p-1} \alpha^{j+1}\left(\Phi_{p j} x^{k}+\phi_{p j}\right) \\
& +\alpha C_{p+1}^{\prime}\left(z_{p+1}-C_{p+1} x^{k}\right)-\alpha C_{p+1}^{\prime} C_{p+1}\left(\psi_{p}-x^{k}\right) \\
= & x^{k}+\alpha \sum_{j=1}^{p+1} C_{j}^{\prime}\left(z_{j}-C_{j} x^{k}\right)+\sum_{j=1}^{p-1} \alpha^{j+1}\left(\Phi_{p j} x^{k}+\phi_{p j}\right) \\
& -\alpha C_{p+1}^{\prime} C_{p+1}\left(\alpha \sum_{j=1}^{p} C_{j}^{\prime}\left(z_{j}-C_{j} x^{k}\right)+\sum_{j=1}^{p-1} \alpha^{j+1}\left(\Phi_{p j} x^{k}+\phi_{p j}\right)\right)
\end{aligned}
$$

which is of the form (5.15), thereby completing the induction.
For $i=m$, Eq. (5.15) yields

$$
x^{k+1}=A(\alpha) x^{k}+b(\alpha)
$$

where

$$
\begin{gathered}
A(\alpha)=I-\alpha \sum_{j=1}^{m} C_{j}^{\prime} C_{j}+\sum_{j=1}^{m-1} \alpha^{j+1} \Phi_{m j} \\
b(\alpha)=-\alpha \sum_{j=1}^{m} C_{j}^{\prime} z_{j}+\sum_{j=1}^{m-1} \alpha^{j+1} \phi_{m j}
\end{gathered}
$$

Let us choose $\alpha$ small enough so that the eigenvalues of $A(\alpha)$ are all strictly within the unit circle: this is possible since $\sum_{j=1}^{m} C_{j}^{\prime} C_{j}{ }^{\prime}$ is assumed positive definite and the last term in Eq. (5.17) involves powers of $\alpha$ that are greater or equal to 2 . Define

$$
x(\alpha)=(I-A(\alpha))^{-1} b(\alpha)
$$

Then $b(\alpha)=(I-A(\alpha)) x(\alpha)$, and by substituting this expression in Eq. (5.16), it can be seen that

$$
x^{k+1}-x(\alpha)=A(\alpha)\left(x^{k}-x(\alpha)\right)
$$

from which

$$
x^{k+1}-x(\alpha)=A(\alpha)^{k}\left(x^{0}-x(\alpha)\right), \quad \forall k
$$

Since all the eigenvalues of $A(\alpha)$ are strictly within the unit circle, we have $A(\alpha)^{k} \rightarrow 0$ (Prop. A. 16 in Appendix A), so $x^{k} \rightarrow x(\alpha)$.

---

To prove that $\lim _{\alpha \rightarrow 0} x(\alpha)=x^{*}$, we first calculate $x^{*}$. We set the gradient of $\|g(x)\|^{2}$ to zero, to obtain

$$
\sum_{i=1}^{m} C_{i}^{\prime}\left(C_{i} x^{*}-z_{i}\right)=0
$$

so that

$$
x^{*}=\left(\sum_{i=1}^{m} C_{i}^{\prime} C_{i}\right)^{-1} \sum_{i=1}^{m} C_{i}^{\prime} z_{i}
$$

Then, we use Eq. (5.19) to write $x(\alpha)=(I / \alpha-A(\alpha) / \alpha)^{-1}(b(\alpha) / \alpha)$, and we see from Eqs. (5.17) and (5.18) that

$$
\lim _{\alpha \rightarrow 0} x(\alpha)=\left(\sum_{i=1}^{m} C_{i}^{\prime} C_{i}\right)^{-1} \sum_{i=1}^{m} C_{i}^{\prime} z_{i}=x^{*}
$$

(b) For $i=m$ and $\alpha=\alpha^{k}$, Eq. (5.15) yields

$$
x^{k+1}=x^{k}+\alpha^{k} \sum_{j=1}^{m} C_{j}^{\prime}\left(z_{j}-C_{j} x^{k}\right)+\left(\alpha^{k}\right)^{2} H^{k}\left(x^{k}-x^{*}\right)+\left(\alpha^{k}\right)^{2} h^{k}
$$

where

$$
\begin{gathered}
H^{k}=\sum_{j=1}^{m-1}\left(\alpha^{k}\right)^{j-1} \Phi_{m j} \\
h^{k}=\sum_{j=1}^{m-1}\left(\alpha^{k}\right)^{j-1}\left(\Phi_{m j} x^{*}+\phi_{m j}\right)
\end{gathered}
$$

Using also the expression (5.22) for $x^{*}$, we can write Eq. (5.23) as

$$
x^{k+1}-x^{*}=\left(I-\alpha^{k} \sum_{j=1}^{m} C_{j}^{\prime} C_{j}+\left(\alpha^{k}\right)^{2} H^{k}\right)\left(x^{k}-x^{*}\right)+\left(\alpha^{k}\right)^{2} h^{k}
$$

Let $k$ be large enough so that a positive number $\xi$ is a lower bound to the minimum eigenvalue of $\sum_{j=1}^{m} C_{j}^{\prime} C_{j}-\alpha^{k} H^{k}$, while $\alpha^{k} \xi<1$. Let also $\delta$ be an upper bound to $\left\|h^{k}\right\|$. Then from Eq. (5.26) we obtain

$$
\left\|x^{k+1}-x^{*}\right\| \leq\left(1-\alpha^{k} \xi\right)\left\|x^{k}-x^{*}\right\|+\left(\alpha^{k}\right)^{2} \delta
$$

Lemma 1.5.1, together with this relation and the assumption $\sum_{k=0}^{\infty} \alpha^{k}=$ $\infty$, imply that $x^{k} \rightarrow x^{*}$. Q.E.D.

In the case where the data blocks are nonlinear, stationarity of the limit points of sequences $\left\{x^{k}\right\}$ generated by the incremental gradient method

---

has been shown under certain assumptions (including Lipschitz continuity of the data block gradients) for the case of the stepsize $\alpha^{k}=\gamma /(k+1)$, where $\gamma$ is a positive scalar [see [MaS94], which in addition analyzes the heavy ball-like method (5.11)]. Contrary to the case of linear data blocks, $\gamma$ may have to be chosen sufficiently small to guarantee boundedness of $\left\{x^{k}\right\}$. The convergence proof is similar to the one of the preceding proposition, but it is technically more involved. In the case of a constant stepsize and nonlinear data blocks, it is also possible to show a result analogous to Prop. 1.5.1(a), but again the proof is technically complex.

# 1.5.3 Incremental Forms of the Gauss-Newton Method* 

We now consider incremental versions of the Gauss-Newton method. An example of such a method starts with some $x^{0}$, then updates $x$ via a Gauss-Newton-like iteration aimed at minimizing

$$
\left\|g_{1}(x)\right\|^{2}
$$

then updates $x$ via a Gauss-Newton-like iteration aimed at minimizing

$$
\lambda\left\|g_{1}(x)\right\|^{2}+\left\|g_{2}(x)\right\|^{2}
$$

where $\lambda$ is a scalar with

$$
0 \leq \lambda \leq 1
$$

and similarly continues, with the $i$ th step consisting of a Gauss-Newton-like iteration aimed at minimizing the weighted partial sum

$$
\sum_{j=1}^{i} \lambda^{i-j}\left\|g_{j}(x)\right\|^{2}
$$

Once the entire data set is processed, the cycle is restarted. The parameter $\lambda$ determines the influence of old data blocks on new estimates. Generally, as $\lambda$ tends towards zero, the effect of old data blocks is discounted faster, and successive estimates produced by the method tend to change more rapidly. Thus one may obtain a faster rate of progress of the method when $\lambda<1$, and this is the main motivation for considering this case.

## The Kalman Filter for Linear Least Squares

When the data blocks are linear functions, it takes a single pure Gauss-Newton iteration to find the least squares estimate, and it turns out that this iteration can be implemented with an incremental algorithm known as the Kalman filter. This algorithm has many important applications in control and communication theory, and has been studied extensively. We first develop this algorithm for linear data blocks, and we then extend it to the case of nonlinear data blocks.

---

Suppose that the functions $g_{i}$ are linear of the form

$$
g_{i}(x)=z_{i}-C_{i} x
$$

where $z_{i} \in \Re^{r_{i}}$ are given vectors and $C_{i}$ are given $r_{i} \times n$ matrices. In other words, we are trying to fit a linear model to the set of measurements $z_{1}, \ldots, z_{m}$. Let us consider the incremental method that sequentially generates the vectors

$$
\psi_{i}=\arg \min _{x \in \Re^{n}} \sum_{j=1}^{i} \lambda^{i-j}\left\|z_{j}-C_{j} x\right\|^{2}, \quad i=1, \ldots, m
$$

Then, for $\lambda=1$, the least squares solution is obtained at the last step as

$$
x^{*}=\psi_{m}
$$

Furthermore, the method can be conveniently implemented as shown by the following proposition:

Proposition 1.5.2: (Kalman Filter) Assuming that the matrix $C_{1}^{\prime} C_{1}$ is positive definite, the least squares estimates

$$
\psi_{i}=\arg \min _{x \in \Re^{n}} \sum_{j=1}^{i} \lambda^{i-j}\left\|z_{j}-C_{j} x\right\|^{2}, \quad i=1, \ldots, m
$$

can be generated by the algorithm

$$
\psi_{i}=\psi_{i-1}+H_{i}^{-1} C_{i}^{\prime}\left(z_{i}-C_{i} \psi_{i-1}\right), \quad i=1, \ldots, m
$$

where $\psi_{0}$ is an arbitrary vector, and the positive definite matrices $H_{i}$ are generated by

$$
H_{i}=\lambda H_{i-1}+C_{i}^{\prime} C_{i}, \quad i=1, \ldots, m
$$

with

$$
H_{0}=0
$$

More generally, for all integers $i$ and $\bar{i}$ with $1 \leq \bar{i}<i \leq m$ we have

$$
\psi_{i}=\psi_{\bar{i}}+H_{i}^{-1} \sum_{j=\bar{i}+1}^{i} \lambda^{i-j} C_{j}^{\prime}\left(z_{j}-C_{j} \psi_{\bar{i}}\right)
$$

---

Proof: We first establish the result for the case of two data blocks in the following lemma:

Lemma 1.5.2: Let $\zeta_{1}, \zeta_{2}$ be given vectors, and $\Gamma_{1}, \Gamma_{2}$ be given matrices such that $\Gamma_{1}^{\prime} \Gamma_{1}$ is positive definite. Then the vectors

$$
\psi_{1}=\arg \min _{x \in \Re^{n}}\left\|\zeta_{1}-\Gamma_{1} x\right\|^{2}
$$

and

$$
\psi_{2}=\arg \min _{x \in \Re^{n}}\left\{\left\|\zeta_{1}-\Gamma_{1} x\right\|^{2}+\left\|\zeta_{2}-\Gamma_{2} x\right\|^{2}\right\}
$$

are also given by

$$
\psi_{1}=\psi_{0}+\left(\Gamma_{1}^{\prime} \Gamma_{1}\right)^{-1} \Gamma_{1}^{\prime}\left(\zeta_{1}-\Gamma_{1} \psi_{0}\right)
$$

and

$$
\psi_{2}=\psi_{1}+\left(\Gamma_{1}^{\prime} \Gamma_{1}+\Gamma_{2}^{\prime} \Gamma_{2}\right)^{-1} \Gamma_{2}^{\prime}\left(\zeta_{2}-\Gamma_{2} \psi_{1}\right)
$$

where $\psi_{0}$ is an arbitrary vector.

Proof: By carrying out the minimization in Eq. (5.34), we obtain

$$
\psi_{1}=\left(\Gamma_{1}^{\prime} \Gamma_{1}\right)^{-1} \Gamma_{1}^{\prime} \zeta_{1}
$$

yielding for any $\psi_{0}$,

$$
\psi_{1}=\psi_{0}-\left(\Gamma_{1}^{\prime} \Gamma_{1}\right)^{-1} \Gamma_{1}^{\prime} \Gamma_{1} \psi_{0}+\left(\Gamma_{1}^{\prime} \Gamma_{1}\right)^{-1} \Gamma_{1}^{\prime} \zeta_{1}
$$

from which the desired Eq. (5.36) follows.
Also, by carrying out the minimization in Eq. (5.35), we obtain

$$
\psi_{2}=\left(\Gamma_{1}^{\prime} \Gamma_{1}+\Gamma_{2}^{\prime} \Gamma_{2}\right)^{-1}\left(\Gamma_{1}^{\prime} \zeta_{1}+\Gamma_{2}^{\prime} \zeta_{2}\right)
$$

or equivalently, using also Eq. (5.38),

$$
\begin{aligned}
\left(\Gamma_{1}^{\prime} \Gamma_{1}+\Gamma_{2}^{\prime} \Gamma_{2}\right) \psi_{2} & =\Gamma_{1}^{\prime} \zeta_{1}+\Gamma_{2}^{\prime} \zeta_{2} \\
& =\Gamma_{1}^{\prime} \Gamma_{1} \psi_{1}+\Gamma_{2}^{\prime} \zeta_{2} \\
& =\left(\Gamma_{1}^{\prime} \Gamma_{1}+\Gamma_{2}^{\prime} \Gamma_{2}\right) \psi_{1}-\Gamma_{2}^{\prime} \Gamma_{2} \psi_{1}+\Gamma_{2}^{\prime} \zeta_{2}
\end{aligned}
$$

from which, by multiplying both sides with $\left(\Gamma_{1}^{\prime} \Gamma_{1}+\Gamma_{2}^{\prime} \Gamma_{2}\right)^{-1}$, the desired Eq. (5.37) follows. Q.E.D.

---

Proof of Prop. 1.5.2: Equation (5.33) follows by applying Lemma 1.5.2 with the correspondences $\psi_{0} \sim \psi_{0}, \psi_{1} \sim \psi_{\bar{i}}, \psi_{2} \sim \psi_{i}$, and

$$
\begin{gathered}
\zeta_{1} \sim\left(\begin{array}{c}
\sqrt{\lambda^{i-1}} z_{1} \\
\vdots \\
\sqrt{\lambda^{i-\bar{i}} z_{\bar{i}}}
\end{array}\right), \quad \Gamma_{1} \sim\left(\begin{array}{c}
\sqrt{\lambda^{i-1}} C_{1} \\
\vdots \\
\sqrt{\lambda^{i-\bar{i}} C_{\bar{i}}}
\end{array}\right) \\
\zeta_{2} \sim\left(\begin{array}{c}
\sqrt{\lambda^{i-\bar{i}-1}} z_{\bar{i}+1} \\
\vdots \\
z_{i}
\end{array}\right), \quad \Gamma_{2} \sim\left(\begin{array}{c}
\sqrt{\lambda^{i-\bar{i}-1}} C_{\bar{i}+1} \\
\vdots \\
C_{i}
\end{array}\right)
\end{gathered}
$$

and by carrying out the straightforward algebra. Equation (5.31) is the special case of Eq. (5.33) corresponding to $\bar{i}=i-1$. Q.E.D.

Note that the positive definiteness assumption on $C_{1}^{\prime} C_{1}$ in Prop. 1.5.2 is needed to guarantee that the first matrix $H_{1}$ is positive definite and hence invertible; then the positive definiteness of the subsequent matrices $H_{2}, \ldots, H_{m}$ follows from Eq. (5.32). As a practical matter, it is possible to guarantee the positive definiteness of $C_{1}^{\prime} C_{1}$ by lumping a sufficient number of measurements into the first data block ( $C_{1}$ should contain $n$ linearly independent columns). An alternative is to redefine $\psi_{i}$ as

$$
\psi_{i}=\arg \min _{x \in \Re^{n}}\left\{\delta \lambda^{i}\left\|x-\psi_{0}\right\|^{2}+\sum_{j=1}^{i} \lambda^{i-j}\left\|z_{j}-C_{j} x\right\|^{2}\right\}, \quad i=1, \ldots, m
$$

where $\delta$ is a small positive scalar. Then it can be seen from the proof of Prop. 1.5.2 that $\psi_{i}$ is generated by the same equations (5.31) and (5.32), except that the initial condition $H_{0}=0$ is replaced by

$$
H_{0}=\delta I
$$

so that $H_{1}=\lambda \delta I+C_{1}^{\prime} C_{1}$ is positive definite even if $C_{1}^{\prime} C_{1}$ is not. Note, however, that in this case, the last estimate $C_{m}$ is only approximately equal to the least squares estimate $x^{*}$, even if $\lambda=1$ (the approximation error depends on the size of $\delta$ ).

# The Extended Kalman Filter 

Consider now the general case where the data blocks $g_{i}$ are nonlinear. Then a generalization of the Kalman filter, known as the extended Kalman filter (EKF for short), can be used. Like the Gauss-Newton method, it involves linearization of the data blocks and solution of linear least squares problems. However, these problems are solved incrementally using the Kalman filtering algorithm, and the linearization of each data block is done

---

at the latest iterafe that is available when this data block is processed. In particular, a cycle through the data set of the algorithm sequentially generates the vectors

$$
v_{i}=\arg \min _{x \in \mathbb{R}^{n}} \sum_{j=1}^{i} \lambda^{i-x}\left\|\tilde{g}_{j}\left(x, v_{j-1}\right)\right\|^{2}, \quad i=1, \ldots, m
$$

where $\tilde{g}_{j}\left(x, \psi_{j-1}\right)$ are the linearized functions

$$
\tilde{g}_{j}\left(x, \psi_{j-1}\right)=g_{j}\left(\psi_{j-1}\right)+\nabla g_{j}\left(\psi_{j-1}\right)^{\prime}\left(x-\psi_{j-1}\right)
$$

and $v_{0}$ is an initial estimate of $x$. Using the formulas (5.31) and (5.32) of Prop. 1.5.2 with the identifications

$$
z_{i}=g_{i}\left(\psi_{i-1}\right)-\nabla g_{i}\left(\psi_{i-1}\right)^{\prime} \psi_{i-1}, \quad C_{i}=-\nabla g_{i}\left(\psi_{i-1}\right)^{\prime}
$$

this algorithm can be written in the incremental form

$$
\psi_{i}=\psi_{i-1}-H_{i}^{-1} \nabla g_{i}\left(\psi_{i-1}\right) g_{i}\left(\psi_{i-1}\right), \quad i=1, \ldots, m
$$

where the matrices $H_{i}$ are generated by

$$
H_{i}=\lambda H_{i-1}+\nabla g_{i}\left(\psi_{i-1}\right) \nabla g_{i}\left(\psi_{i-1}\right)^{\prime}, \quad i=1, \ldots, m
$$

with

$$
H_{0}=0
$$

To contrast the EKF with the pure form of the Gauss-Newton method (unit stepsize), note that a single iteration of the latter can be written as

$$
x^{k+1}=\arg \min _{x \in \mathbb{R}^{m}} \sum_{i=1}^{m}\left\|\tilde{g}_{i}\left(x, x^{k}\right)\right\|^{2}
$$

Using the formulas of Prop. 1.5.2 with the identifications

$$
z_{i}=g_{i}\left(x^{k}\right)-\nabla g_{i}\left(x^{k}\right)^{\prime} x^{k}, \quad C_{i}=-\nabla g_{i}\left(x^{k}\right)^{\prime}
$$

we can generate $x^{k+1}$ by an incremental algorithm as

$$
x^{k+1}=\bar{\psi}_{m}
$$

where
$\bar{\psi}_{i}=\bar{\psi}_{i-1}-\bar{H}_{i}^{-1} \nabla g_{i}\left(x^{k}\right)\left(g_{i}\left(x^{k}\right)+\nabla g_{i}\left(x^{k}\right)^{\prime}\left(\bar{\psi}_{i-1}-x^{k}\right)\right), \quad i=1, \ldots, m$,
$\bar{\psi}_{0}=x^{k}$, and the matrices $\bar{H}_{i}$ are generated by

$$
\bar{H}_{i}=\bar{H}_{i-1}+\nabla g_{i}\left(x^{k}\right) \nabla g_{i}\left(x^{k}\right)^{\prime}, \quad i=1, \ldots, m
$$

with

$$
\bar{H}_{0}=0
$$

Thus, by comparing Eqs. (5.41)-(5.43) with Eqs. (5.45)-(5.47), we see that, if $\lambda=1$, a cycle of the EKF through the data set differs from a pure Gauss-Newton iteration only in that the linearization of the data blocks $g_{i}$ is done at the corresponding current estimates $\psi_{i-1}$ rather than at the estimate $x^{k}$ available at the start of the cycle.

---

# Convergence Issues for the Extended Kalman Filter 

We have considered so far a single cycle of the EKF. To obtain an algorithm that cycles through the data set multiple times, we can simply create a larger data set by concatenating multiple copies of the original data set, that is, by forming what we refer to as the extended data set

$$
\left(g_{1}, g_{2}, \ldots, g_{m}, g_{1}, g_{2}, \ldots, g_{m}, g_{1}, g_{2}, \ldots\right)
$$

The EKF when applied to the extended data set often works well in practice, although its convergence properties as the number of cycles increases to infinity have not been fully investigated, particularly when $\lambda<1$.

The basic reason why the algorithm works is that asymptotically it resembles a gradient method with diminishing stepsize of the type described in Section 1.2. To get a sense of this, assume that the EKF is applied to the extended data set (5.48) with $\lambda=1$. Let us denote by $x^{k}$ the iterate at the end of the $k$ th cycle through the data set, that is,

$$
x^{k}=\psi_{k m}, \quad k=1,2, \ldots
$$

Then by using Eq. (5.33) with $i=(k+1) m$ and $\bar{i}=k m$, we obtain

$$
x^{k+1}=x^{k}-H_{(k+1) m}^{-1}\left(\sum_{i=1}^{m} \nabla g_{i}\left(\psi_{k m+i-1}\right) g_{i}\left(\psi_{k m+i-1}\right)\right)
$$

Now $H_{(k+1) m}$ grows roughly in proportion to $k+1$ because, by Eq. (5.32), we have

$$
H_{(k+1) m}=\sum_{j=0}^{k} \sum_{i=1}^{m} \nabla g_{i}\left(\psi_{j m+i-1}\right) \nabla g_{i}\left(\psi_{j m+i-1}\right)^{\prime}
$$

It is therefore reasonable to expect that the method tends to make slow progress when $k$ is large, which means that the vectors $\psi_{k m+i-1}$ in Eq. (5.49) are roughly equal to $x^{k}$. Thus for large $k$, the sum in the righthand side of Eq. (5.49) is roughly equal to the gradient $\nabla g\left(x^{k}\right) g\left(x^{k}\right)$, while from Eq. (5.50), $H_{(k+1) m}$ is roughly equal to $(k+1) \nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}$, where $g=\left(g_{1}, g_{2}, \ldots, g_{m}\right)$ is the original data set. It follows that for large $k$, the EKF iteration (5.49) can be written approximately as

$$
x^{k+1} \approx x^{k}-\frac{1}{k+1}\left(\nabla g\left(x^{k}\right) \nabla g\left(x^{k}\right)^{\prime}\right)^{-1} \nabla g\left(x^{k}\right) g\left(x^{k}\right)
$$

that is, as an approximate Gauss-Newton iteration with diminishing stepsize.

When $\lambda<1$, the matrix $H_{i}^{-1}$ generated by the EKF recursion (5.42) will typically not diminish to zero, and $\left\{x^{k}\right\}$ may not converge to a stationary point of $\sum_{i=1}^{m} \lambda^{m-i}\left\|g_{i}(x)\right\|^{2}$. Furthermore, as the following example shows, the sequences $\left\{\psi_{k m+i}\right\}$ produced by the EKF using Eq. (5.41), may converge to different limits for different $i$ :

---

# Example 5.6 

Consider the case where there are two data blocks, $g_{1}(x)=x-c_{1}$ and $g_{2}(x)=$ $x-c_{2}$, where $c_{1}$ and $c_{2}$ are given scalars. Each cycle of the EKF consists of two steps. At the second step of the $k$ th cycle, we minimize

$$
\sum_{i=1}^{k}\left(\lambda^{2 i-1}\left(x-c_{1}\right)^{2}+\lambda^{2 i-2}\left(x-c_{2}\right)^{2}\right)
$$

which is equal to the following scalar multiple of $\lambda\left(x-c_{1}\right)^{2}+\left(x-c_{2}\right)^{2}$,

$$
\left(1+\lambda^{2}+\cdots+\lambda^{2 k-2}\right)\left(\lambda\left(x-c_{1}\right)^{2}+\left(x-c_{2}\right)^{2}\right)
$$

Thus at the second step, we obtain the minimizer of $\lambda\left(x-c_{1}\right)^{2}+\left(x-c_{2}\right)^{2}$,

$$
\psi_{2 k}=\frac{\lambda c_{1}+c_{2}}{\lambda+1}
$$

At the first step of the $k$ th cycle, we minimize

$$
\left(x-c_{1}\right)^{2}+\lambda \sum_{i=1}^{k-1}\left(\lambda^{2 i-1}\left(x-c_{1}\right)^{2}+\lambda^{2 i-2}\left(x-c_{2}\right)^{2}\right)
$$

which is equal to the following scalar multiple of $\left(x-c_{1}\right)^{2}+\lambda\left(x-c_{2}\right)^{2}$

$$
\left(1+\lambda^{2}+\cdots+\lambda^{2 k-4}\right)\left(\left(x-c_{1}\right)^{2}+\lambda\left(x-c_{2}\right)^{2}\right)
$$

plus the diminishing term $\lambda^{2 k-2}\left(x-c_{1}\right)^{2}$. Thus at the first step, we obtain approximately (for large $k$ ) the minimizer of $\left(x-c_{1}\right)^{2}+\lambda\left(x-c_{2}\right)^{2}$,

$$
\psi_{2 k-1} \approx \frac{c_{1}+\lambda c_{2}}{1+\lambda}
$$

We see therefore that within each cycle, there is an oscillation around the minimizer $\left(c_{1}+c_{2}\right) / 2$ of $\left(x-c_{1}\right)^{2}+\left(x-c_{2}\right)^{2}$. The size of the oscillation diminishes as $\lambda$ approaches 1 .

Generally, for a nonlinear least squares problem, the convergence rate tends to be faster when $\lambda<1$ than when $\lambda=1$, essentially because the implicit stepsize does not diminish as in the case $\lambda=1$. For this reason, a hybrid method that uses a different value of $\lambda$ within each cycle may work best in practice. One may start with a relatively small $\lambda$ to attain a fast initial rate of convergence, and then progressively increase $\lambda$ towards 1 in order to attain high solution accuracy. An example that the reader may try is given in Exercise 5.7. The following proposition shows convergence for the case where $\lambda$ tends to 1 at a sufficiently fast rate.

---

Proposition 1.5.3: Assume that $\nabla g_{i}(x)$ has full rank for all $x$ and $i=1, \ldots, m$, and that for some $L>0$, we have

$$
\left\|\nabla g_{i}(x) g_{i}(x)-\nabla g_{i}(y) g_{i}(y)\right\| \leq L\|x-y\|, \quad \forall x, y \in \Re^{n}, i=1, \ldots, m
$$

Assume also that there is a constant $c>0$ such that the scalar $\lambda$ used in the updating formula (5.42) within the $k$ th cycle, call it $\lambda(k)$, satisfies

$$
0 \leq 1-(\lambda(k))^{m} \leq \frac{c}{k}, \quad \forall k=1,2, \ldots
$$

Then if the EKF applied to the extended data set (5.48) generates a bounded sequence of vectors $\psi_{i}$, and each of the limit points of $\left\{x^{k}\right\}$ is a stationary point of the least squares problem.

Proof: (Abbreviated; for a complete proof see [Ber94]) We have using the Kalman filter recursion (5.33) that $x^{k}$ satisfies

$$
x^{k+1}=x^{k}-H_{(k+1) m}^{-1}\left(\sum_{i=1}^{m}(\lambda(k))^{m-i} \nabla g_{i}\left(\psi_{k m+i-1}\right) g_{i}\left(\psi_{k m+i-1}\right)\right)
$$

It can be shown using the rank assumption on $\nabla g_{i}(x)$, the growth assumption on $\lambda^{k}$, the boundedness of $\left\{x^{k}\right\}$, and the preceding analysis that the eigenvalues of the matrices $H_{k m}$ are within an interval $\left[c_{1} k, c_{2} k\right]$, where $c_{1}$ and $c_{2}$ are some positive constants (see Exercise 5.9). The proof then follows the line of argument of the convergence proof of gradient methods with diminishing stepsize (Prop. 1.2.4 in Section 1.2). Q.E.D.

# E XERCISES 

## 5.1

Consider the least squares problem (5.1) for the case where $m<n$.
(a) Show that the Hessian matrix at any optimal solution is singular.
(b) Consider the case where $g$ is linear and of the form $g(x)=z-A x$, where $A$ is an $m \times n$ matrix. Show that there are infinitely many optimal solutions. Show also that if $A$ has linearly independent rows, $x^{*}=A^{\prime}\left(A A^{\prime}\right)^{-1} z$ is one of these solutions.

---

# 5.2 [Luo91] 

Consider the least squares problem

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2}\left\{\left(z_{1}-x\right)^{2}+\left(z_{2}-x\right)^{2}\right\} \\
& \text { subject to } x \in \Re
\end{aligned}
$$

and the incremental gradient algorithm that generates $x^{k+1}$ from $x^{k}$ according to

$$
x^{k+1}=y^{k}-\alpha\left(y^{k}-z_{2}\right)
$$

where

$$
y^{k}=x^{k}-\alpha\left(x^{k}-z_{1}\right)
$$

and $\alpha$ is positive stepsize. Assuming that $\alpha<1$, show that $\left\{x^{k}\right\}$ and $\left\{y^{k}\right\}$ converge to limits $x(\alpha)$ and $y(\alpha)$, respectively. However, unless $z_{1} \neq z_{2}$, $x(\alpha)$ and $y(\alpha)$ are neither equal to each other, nor equal to the least squares solution $x^{*}=\left(z_{1}+z_{2}\right) / 2$. Consistently with Prop. 1.5.1, verify that

$$
\lim _{\alpha \rightarrow 0} x(\alpha)=\lim _{\alpha \rightarrow 0} y(\alpha)=x^{*}
$$

### 5.3 (Computational Problem)

Consider the least squares cost function of Fig. 1.5.1 for the case where the five data pairs $\left(y_{s}, z_{s}\right)$ are $(1.165,1),(0.626,-1),(0.075,-1),(0.351,1)$, $(-0.696,1)$ (these correspond to the three-dimensional plot of Fig. 1.5.1). Minimize this cost function using appropriate versions of the Gauss-Newton method, the Extended Kalman Filter, the incremental gradient method, and its heavy ball version (5.11).

## 5.4 (A Generalized Incremental Gradient Method [Ber95b])

The purpose of this exercise is to embed the incremental gradient method and the method of steepest descent within a one-parameter family of methods for the least squares problem. In contrast with the incremental gradient method, some of the methods in this family have a generically linear convergence rate. For a fixed $\mu \geq 0$, define

$$
\xi_{i}(\mu)=\frac{1}{1+\mu+\cdots+\mu^{m-i}}, \quad i=1, \ldots, m
$$

Consider a method which given $x^{k}$, generates $x^{k+1}$ according to $x^{k+1}=\psi_{m}$, where

$$
\begin{gathered}
\psi_{i}=x^{k}-\alpha^{k} h_{i} \\
h_{i}=\mu h_{i-1}+\sum_{j=1}^{i} \xi_{j}(\mu) \nabla g_{j}\left(\psi_{j-1}\right) g_{j}\left(\psi_{j-1}\right)
\end{gathered}
$$

---

with the initial condition $\dot{h}_{0}=0$.
(a) Show that when $\mu=0$, the method coincides with the incremental gradient method, and when $\mu \rightarrow \infty$, the method approaches the steepest descent method.
(b) Show that parts (a) and (b) of Prop. 1.5.1, which were proved for the case $\mu=0$, hold as stated when $\mu>0$ as well.
(c) Suppose that in the $k$ th iteration, a $k$-dependent value of $\mu$, say $\mu(k)$, and a constant stepsize $\alpha^{k}=\alpha$ are used. Under the assumptions of Prop. 1.5.1, show that if for some $q>1$ and all $k$ greater than some index $\bar{k}$, we have $\mu(k) \geq q^{k}$, then there exists $\bar{\alpha}>0$ such that for all $\alpha \in(0, \bar{\alpha}]$, the iteration converges linearly to the optimal solution $x^{*}$.

# 5.5 (Gradient Method for Infinitely Many Data Blocks) 

Consider the gradient method

$$
x^{k+1}=x^{k}-\alpha \nabla f_{k}\left(x^{k}\right) \quad k=0,1, \ldots
$$

where $f_{0}, f_{1}, \ldots$, are quadratic functions with eigenvalues lying within some interval $[\gamma, \Gamma]$, where $\gamma>0$. Suppose that for a given $\epsilon>0$, there is a vector $x^{*}$ such that

$$
\left\|\nabla f_{k}\left(x^{*}\right)\right\| \leq \epsilon, \quad \forall k=0,1, \ldots
$$

Show that for all $\alpha$ with $0<\alpha \leq 2 /(\gamma+\Gamma)$, we have

$$
\limsup _{k \rightarrow \infty}\left\|x^{k}-x^{*}\right\| \leq \frac{2 \epsilon}{\gamma}
$$

Hint: Let $Q_{k}$ be the positive definite symmetric matrix corresponding to $f_{k}$, and write

$$
x^{k+1}-x^{*}=\left(I-\alpha Q_{k}\right)\left(x^{k}-x^{*}\right)-\alpha \nabla f_{k}\left(x^{*}\right)
$$

Use this relation to show that

$$
\left\|x^{k}-x^{*}\right\|>\frac{2 \epsilon}{\gamma} \quad \Rightarrow \quad\left\|x^{k+1}-x^{*}\right\|<\left(1-\frac{\alpha \gamma}{2}\right)\left\|x^{k}-x^{*}\right\|
$$

while

$$
\left\|x^{k}-x^{*}\right\| \leq \frac{2 \epsilon}{\gamma} \quad \Rightarrow \quad\left\|x^{k+1}-x^{*}\right\| \leq \frac{2 \epsilon}{\gamma}
$$

5.6

Show that if the sequence $\left\{x^{k}\right\}$ produced by the EKF with $\lambda=1$ converges to some $\bar{x}$, then $\bar{x}$ must be a stationary point of $\sum_{i=1}^{m}\left\|g_{i}(x)\right\|^{2}$.

---

# 5.7 

Consider the EKF of Eqs. (5.41)-(5.43) for the data set $g_{1}(x)=x^{2}-2$ and $g_{2}(x)=x^{2}$, where $x \in \Re$.
(a) Verify that each cycle of the EKF consists of executing sequentially the iterations

$$
\begin{gathered}
H_{y}:=\lambda H_{x}+(2 x)^{2} \\
y:=x-H_{y}^{-1}(2 x)\left(x^{2}-2\right) \\
H_{x}:=\lambda H_{y}+(2 y)^{2} \\
x:=y-H_{x}^{-1}(2 y)\left(y^{2}\right)
\end{gathered}
$$

(b) Write a computer program that implements the EKF, and compare the limit of $(x+y) / 2$ with the least squares solutions, which are 1 and -1 , for different values of $\lambda$. Try also a version of the method where $\lambda$ is progressively increased towards 1.

### 5.8 (Incremental Gauss-Newton Method with Restart)

Consider a version of the EKF where the matrix $H$ is reset to 0 at the end of each cycle. In particular, this method, given $x^{k}$, generates $x^{k+1}$ according to $x^{k+1}=\psi_{m}$, where

$$
\begin{aligned}
\psi_{i}=\psi_{i-1}-\alpha^{k} H_{i}^{-1} \nabla g_{i}\left(\psi_{i-1}\right) g_{i}\left(\psi_{i-1}\right), & i=1, \ldots, m \\
H_{i}=H_{i-1}+\nabla g_{i}\left(\psi_{i-1}\right) \nabla g_{i}\left(\psi_{i-1}\right)^{\prime}, & i=1, \ldots, m
\end{aligned}
$$

with the initial conditions

$$
\psi_{0}=x^{k}, \quad H_{0}=0
$$

(a) Show that if the data blocks $g_{i}$ are linear and $\alpha^{k}$ is for all $k$ equal to a constant $\alpha$ with $0<\alpha<2$, then $\left\{x^{k}\right\}$ converges to the optimal solution.
(b) Show that if $\sum_{k=0}^{\infty} \alpha^{k}=\infty, \sum_{k=0}^{\infty}\left(\alpha^{k}\right)^{2}<\infty$, and the generated sequence $\left\{x^{k}\right\}$ is bounded, then every limit point of $\left\{x^{k}\right\}$ is a stationary point of the least squares problem.

## 5.9

Under the assumptions of Prop. 1.5.3, show that the eigenvalues of the matrices $H_{k m}$ lie within an interval $\left[c_{1} k, c_{2} k\right]$, where $c_{1}$ and $c_{2}$ are some positive constants. Hint: Let $X$ be a compact set containing all vectors $\psi_{i}$ generated by the algorithm, and let $B$ and $b$ be an upper bound and a lower bound, respectively, on the eigenvalues of $\nabla g_{i}(x) \nabla g_{i}(x)^{\prime}$ as $x$ ranges over $X$. Show that all eigenvalues of $H_{k m}$ are less or equal to $k m B$. If $v_{k}$ is the smallest eigenvalue of $H_{k m}$, show that $v_{k+1} \geq(1-c / k) v_{k}+m(1-c / k) b$, and use this relation to show that $v_{k} \geq k \gamma$ for a sufficiently small but positive value of $\gamma$.

---

# 1.6 CONJUGATE DIRECTION METHODS 

Conjugate direction methods are motivated by a desire to accelerate the convergence rate of steepest descent, while avoiding the overhead associated with Newton's method. They were originally developed for solving the quadratic problem

$$
\begin{array}{ll}
\operatorname{minimize} & f(x)=\frac{1}{2} x^{\prime} Q x-b^{\prime} x \\
\text { subject to } & x \in \mathbb{R}^{n}
\end{array}
$$

where $Q$ is positive definite, or equivalently, for solving the linear system

$$
Q x=b
$$

They can also be used for solution of the more general system $A x=b$, where $A$ is invertible but not positive definite, after conversion to the positive definite system $A^{\prime} A x=A^{\prime} b$.

Conjugate direction methods can solve these problems after at most $n$ iterations but they are best viewed as iterative methods, since usually fewer than $n$ iterations are required to attain a sufficiently accurate solution, particularly when $n$ is large. They can also be used to solve nonquadratic optimization problems. For such problems, they do not in general terminate after a finite number of iterations, but still, when properly implemented, they have attractive convergence and rate of convergence properties. We will first develop the methods for quadratic problems and then discuss their application to more general problems.

Given a positive definite $n \times n$ matrix $Q$, we say that a set of nonzero vectors $d^{1}, \ldots, d^{k}$ are $Q$-conjugate, if

$$
d^{i^{\prime}} Q d^{j}=0, \quad \text { for all } i \text { and } j \text { such that } i \neq j
$$

If $d^{1}, \ldots, d^{k}$ are $Q$-conjugate, then they are linearly independent, since if one of these vectors, say $d^{k}$, were expressed as a linear combination of the others,

$$
d^{k}=\alpha^{1} d^{1}+\cdots+\alpha^{k-1} d^{k-1}
$$

then by multiplication with $d^{k^{\prime}} Q$ we would obtain using the $Q$-conjugacy of $d^{k}$ and $d^{j}, j=1, \ldots, k-1$,

$$
d^{k^{\prime}} Q d^{k}=\alpha^{1} d^{k^{\prime}} Q d^{1}+\cdots+\alpha^{k-1} d^{k^{\prime}} Q d^{k-1}=0
$$

which is impossible since $d^{k} \neq 0$ and $Q$ is positive definite.
For a given set of $n Q$-conjugate directions $d^{0}, \ldots, d^{n-1}$, the corresponding conjugate direction method for unconstrained minimization of the quadratic function

$$
f(x)=\frac{1}{2} x^{\prime} Q x-b^{\prime} x
$$

---

is given by

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}, \quad k=0, \ldots, n-1
$$

where $x^{0}$ is an arbitrary starting vector and $\alpha^{k}$ is obtained by the line minimization rule

$$
f\left(x^{k}+\alpha^{k} d^{k}\right)=\min _{\alpha} f\left(x^{k}+\alpha d^{k}\right)
$$

The principal result about conjugate direction methods is that successive iterates minumize $f$ over a progressively expanding linear manifold that eventually includes the global minimum of $f$. In particular, for each $k$, $x^{k+1}$ minimizes $f$ over the linear manifold passing through $x^{0}$ and spanned by the conjugate directions $d^{0}, \ldots, d^{k}$, that is,

$$
x^{k+1}=\arg \min _{x \in M^{k}} f(x)
$$

where

$$
\begin{aligned}
M^{k} & =\left\{x \mid x=x^{0}+v, v \in\left(\text { subspace spanned by } d^{0}, \ldots, d^{k}\right)\right\} \\
& =x^{0}+\left(\text { subspace spanned by } d^{0}, \ldots, d^{k}\right)
\end{aligned}
$$

In particular, $x^{n}$ minimizes $f$ over $\Re^{n}$.
To show this, note that by Eq. (6.6), we have for all $i$

$$
\left.\frac{\partial f\left(x^{i}+\alpha d^{i}\right)}{\partial \alpha}\right|_{\alpha=\alpha^{\prime}}=\nabla f\left(x^{i+1}\right)^{\prime} d^{i}=0
$$

and, for $i=0, \ldots, k-1$,

$$
\begin{aligned}
\nabla f\left(x^{k+1}\right)^{\prime} d^{i} & =\left(Q x^{k+1}-b\right)^{\prime} d^{i} \\
& =\left(x^{i+1}+\sum_{j=i+1}^{k} \alpha^{j} d^{j}\right)^{\prime} Q d^{i}-b^{\prime} d^{i} \\
& =x^{i+1^{\prime}} Q d^{i}-b^{\prime} d^{i} \\
& =\nabla f\left(x^{i+1}\right)^{\prime} d^{i}
\end{aligned}
$$

where we have used the conjugacy of $d^{i}$ and $d^{j}, j=i+1, \ldots k$. Combining the last two equations we obtain

$$
\nabla f\left(x^{k+1}\right)^{\prime} d^{i}=0, \quad i=0, \ldots, k
$$

so that

$$
\left.\frac{\partial f\left(x^{0}+\gamma^{0} d^{0}+\cdots+\gamma^{k} d^{k}\right)}{\partial \gamma^{i}}\right|_{\substack{\gamma^{j}=\alpha^{j} \\ j=1 \ldots, k}}=0, \quad i=0, \ldots, k
$$

---

which verifies Eq. (6.7).
It is easy to visualize the expanding manifold minimization property of Eq. (6.7) when $b=0$ and $Q=I$ (the identity matrix). In this case, the equal cost surfaces of $f$ are concentric spheres, and the notion of $Q$ conjugacy reduces to usual orthogonality. By a simple algebraic argument, we see that minimization along $n$ orthogonal directions yields the global minimum of $f$, that is, the center of the spheres. (This becomes evident once we rotate the coordinate system so that the given $n$ orthogonal directions coincide with the coordinate directions.) The case of a general positive definite $Q$ can be reduced to the case where $Q=I$ by means of a scaling transformation. By setting $y=Q^{1 / 2} x$, minimizing $\frac{1}{2} x^{\prime} Q x$ is equivalent to minimizing $\frac{1}{2}\|y\|^{2}$. If $w^{0}, \ldots, w^{n-1}$ are any set of orthogonal nonzero vectors in $\Re^{n}$, the algorithm

$$
y^{k+1}=y^{k}+\alpha^{k} w^{k}, \quad k=0, \ldots, n-1
$$

where

$$
\alpha^{k}=\arg \min _{n} \frac{1}{2}\left\|y^{k}+\alpha w^{k}\right\|^{2}
$$

terminates in at most $n$ steps with $y^{n}=0$. To pass back to the $x$-coordinate system, we multiply Eq. (6.10) by $Q^{-1 / 2}$ and obtain

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}, \quad k=0, \ldots, n-1
$$

where $d^{k}=Q^{-1 / 2} w^{k}$. The orthogonality property of $w^{0}, \ldots, w^{n-1}$, that is, $w^{i^{\prime}} w^{j}=0$ for $i \neq j$, is equivalent to the requirement that the directions $d^{0}, \ldots, d^{n-1}$ be $Q$-conjugate, that is, $d^{i^{\prime}} Q d^{j}=0$ for $i \neq j$.

Thus, using the transformation $y=Q^{1 / 2} x$, we can think of any conjugate direction method for minimizing $\frac{1}{2} x^{\prime} Q x$ as a method that minimizes $\frac{1}{2}\|y\|^{2}$ by successive minimization along $n$ orthogonal directions (see Fig. 1.6.1).

# Generating $Q$-Conjugate Directions 

Given any set of linearly independent vectors $\xi^{0}, \ldots, \xi^{k}$, we can construct a set of mutually $Q$-conjugate directions $d^{0}, \ldots d^{k}$ such that for all $i=0, \ldots, k$, we have
(subspace spanned by $d^{0}, \ldots, d^{i}$ ) $=$ (subspace spanned by $\left.\xi^{0}, \ldots, \xi^{i}\right)$,
using the so called Gmm-Schmidt procedure. Indeed, let us do this recursively, starting with

$$
d^{0}=\xi^{0}
$$

Suppose that, for some $i<k$, we have selected $Q$-conjugate $d^{0}, \ldots, d^{i}$ so that the above property holds. We then take $d^{i+1}$ to be of the form

$$
d^{i+1}=\xi^{i+1}+\sum_{m=0}^{i} c^{(i+1) m} d^{m}
$$

---

![[chapter_21_p121_img26.jpeg]]

Figure 1.6.1. Geometric interpretation of conjugate direction methods in terms of successive minimization along $n$ orthogonal directions. In (a) the function $\|y\|^{2}$ is minimized successively along the directions $w^{0}, \ldots, w^{n-1}$, which are orthogonal in the usual sense ( $w^{i^{\prime}} w^{j}=0$ for $i \neq j$ ). When this process is viewed in the coordinate system of variables $x=Q^{-1 / 2} y$, it yields the conjugate direction method that uses the $Q$-conjugate directions $d^{0}, \ldots, d^{n-1}$ with $d^{i}=Q^{-1 / 2} w^{i}$, as shown in (b).
and choose the coefficients $c^{(i+1) m}$ so that $d^{i+1}$ is $Q$-conjugate to $d^{0}, \ldots, d^{i}$. This will be so if for each $j=0, \ldots, i$,

$$
d^{i+1^{\prime}} Q d^{j}=\xi^{i+1^{\prime}} Q d^{j}+\left(\sum_{m=0}^{i} c^{(i+1) m} d^{m}\right)^{\prime} Q d^{j}=0
$$

Since $d^{0}, \ldots, d^{i}$ are $Q$-conjugate, we have $d^{m^{\prime}} Q d^{j}=0$ if $m \neq j$, and Eq. (6.14) yields

$$
c^{(i+1) j}=-\frac{\xi^{i+1^{\prime}} Q d^{j}}{d^{j^{\prime}} Q d^{j}}, \quad j=0, \ldots, i
$$

Note that the denominator $d^{j^{\prime}} Q d^{j}$ in the above equation is nonzero, since $d^{0}, \ldots d^{i}$ are assumed $Q$-conjugate and are therefore nonzero. Note also that $d^{i+1} \neq 0$, since otherwise from Eqs. (6.11) and (6.13), $\xi^{i+1}$ would be a linear combination of $\xi^{0}, \ldots, \xi^{i}$, contradicting the linear independence of $\xi^{0}, \ldots, \xi^{k}$. Finally, note from Eq. (6.13) that $\xi^{i+1}$ lies in the subspace spanned by $d^{0}, \ldots d^{i+1}$, while $d^{i+1}$ lies in the subspace spanned by $\xi^{0}, \ldots, \xi^{i+1}$, since $d^{0}, \ldots, d^{i}$ and $\xi^{0}, \ldots, \xi^{i}$ span the same space [cf. Eq. (6.11)]. Thus, Eq. (6.11) is satisfied when $i$ is increased to $i+1$ and the Gram-Schmidt procedure defined by Eqs. (6.12), (6.13), and (6.15), has the property claimed. Figure 1.6.2 illustrates the procedure.

It is also worth noting what will happen if the vectors $\xi^{0}, \ldots, \xi^{i}$ are linearly independent, but the next vector $\xi^{i+1}$ is linearly dependent on these vectors. In this case it can be seen (compare also with Fig. 1.6.2) that the equations (6.13) and (6.15) are still valid, but the new vector $d^{i+1}$

---

![[chapter_21_p122_img27.jpeg]]

Figure 1.6.2. Illustration of the Gram-Schmidt procedure for generating $Q$ conjugate directions $d^{0}, \ldots, d^{k}$ from a set of linearly independent vectors $\xi^{0}, \ldots, \xi^{k}$, so that
(subspace spanned by $d^{0}, \ldots, d^{k}$ ) = (subspace spanned by $\xi^{0}, \ldots, \xi^{k}$ ).
Given $d^{0}, \ldots, d^{i-1}$, the $i$ th direction is obtained as $d^{i}=\xi^{i}-\xi^{i}$, where $\xi^{i}$ is a vector on the subspace spanned by $d^{0}, \ldots, d^{i-1}$ (or $\xi^{0}, \ldots, \xi^{i-1}$ ) chosen so that $d^{i}$ is $Q$-conjugate to $d^{0}, \ldots, d^{i-1}$. (It can be shown that the vector $\xi^{i}$ is the projection of $\xi^{i}$ on this subspace with respect to the norm $\|x\|_{Q}=\sqrt{x^{\prime} Q x}$, that is, minimizes $\left\|\xi^{i}-x\right\|_{Q}$ over all $x$ in this subspace; see Exercise 5.1.)
as given by Eq. (6.13) will be zero? We can use this property to construct a set of $Q$-conjugate directions that span the same space as a set of vectors $\xi^{0}, \ldots, \xi^{k}$ that are not a priori known to be linearly independent. This construction can be accomplished with an extended version of the GramSchmidt procedure that generates directions via Eqs. (6.13) and (6.15), but each time the new direction $d^{i+1}$ as given by Eq. (6.13) turns out to be zero, it is simply discarded rather than added to the set of preceding directions.

# The Conjugate Gradient Method 

The most important conjugate direction method, the conjugate gradient method, is obtained by applying the Gram-Schmidt procedure to the gradient vectors $\xi^{0}=-g^{0}, \ldots, \xi^{n-1}=-g^{n-1}$, where we use the notation

$$
g^{k}=\nabla f\left(x^{k}\right)=Q x^{k}-b
$$

Thus the conjugate gradient method is defined by

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}
$$

[^0]
[^0]:    ${ }^{5}$ From Eq. (6.13) and the linear independence of $\xi^{0}, \ldots, \xi^{i}$, it is seen that $d^{i+1}$ can be uniquely expressed as $d^{i+1}=\sum_{m=1}^{i} \gamma^{m} d^{m}$, where $\gamma^{m}$ are some scalars. By multiplying this equation with $d^{2} Q$ and by using the $Q$-conjugacy of the directions $d^{0}, \ldots, d^{i}$ and Eq. (6.14), we see that $\gamma^{m}=0$ for all $m=0, \ldots, i$.

---

where the stepsize $\alpha^{k}$ is obtained by line minimization, and the direction $d^{k}$ is obtained by applying the $k$ th step of the Gram-Schmidt procedure to the vector $-g^{k}$ and the preceding directions $d^{0}, \ldots, d^{k-1}$. In particular, from the Gram-Schmidt equations (6.13) and (6.15), we have

$$
d^{k}=-g^{k}+\sum_{j=0}^{k-1} \frac{g^{k^{\prime}} Q d^{j}}{d^{j^{\prime}} Q d^{j}} d^{j}
$$

Note here that

$$
d^{0}=-g^{0}
$$

and that the method terminates with an optimal solution if $g^{k}=0$. The method also effectively stops if $d^{k}=0$, but we will show that this can only happen if $g^{k}=0$.

The key property of the conjugate gradient method is that the direction formula (6.18) can be greatly simplified. In particular, all but one of the coefficients in the sum of Eq. (6.18) turn out to be zero because, in view of the expanding manifold minimization property, the gradient $g^{k}$ is orthogonal to the subspace spanned by $d^{0}, \ldots d^{k-1}$ [cf. Eq. (6.9)]. We have the following proposition.

Proposition 1.6.1: The directions of the conjugate gradient method are generated by

$$
\begin{gathered}
d^{0}=-g^{0} \\
d^{k}=-g^{k}+\beta^{k} d^{k-1}, \quad k=1, \ldots, n-1
\end{gathered}
$$

where $\beta^{k}$ is given by

$$
\beta^{k}=\frac{g^{k^{\prime}} g^{k}}{g^{k-1^{\prime}} g^{k-1}}
$$

Furthermore, the method terminates with an optimal solution after at most $n$ steps.

Proof: We first use induction to show that all the gradients $g^{k}$ generated up to termination are linearly independent. The result is clearly true for $k=1$. Suppose that the method has not terminated after $k$ steps, and that $g^{0}, \ldots, g^{k-1}$ are linearly independent. Then, since the method is by definition a conjugate direction method, we have
(subspace spanned by $d^{0}, \ldots, d^{k-1}$ ) = (subspace spanned by $g^{0}, \ldots, g^{k-1}$ )
[cf. Eq. (6.11)]. There are two possibilities:
(a) $g^{k}=0$, in which case the method terminates.

---

(b) $g^{k} \neq 0$, in which case the expanding manifold minimization property of the conjugate direction method [cf. Eq. (6.9)] implies that

$$
g^{k} \text { is orthogonal to } d^{0}, \ldots, d^{k-1}
$$

Since the subspaces spanned by $\left(d^{0}, \ldots, d^{k-1}\right)$ and by $\left(g^{0}, \ldots, g^{k-1}\right)$ are the same [cf. Eq. (6.19)], we see that

$$
g^{k} \text { is orthogonal to } g^{0}, \ldots, g^{k-1}
$$

Therefore, $g^{k}$ is linearly independent of $g^{0} \ldots, g^{k-1}$, thus completing the induction.

Since at most $n$ linearly independent gradients can be generated, it follows that the gradient will be zero after at most $n$ iterations and the method will terminate with the minimum of $f$.

To conclude the proof, we use the orthogonality properties (6.20) and (6.21) to verify that the calculation of the coefficients multiplying $d^{j}$ in the Gram-Schmidt formula (6.18) can be simplified as stated in the proposition. We have for all $j$ such that $g^{j} \neq 0$,

$$
g^{j+1}-g^{j}=Q\left(x^{j+1}-x^{j}\right)=\alpha^{j} Q d^{j}
$$

[cf. Eqs. (6.16) and (6.17)]. We note that $\alpha^{j} \neq 0$, since if $\alpha^{j}=0$ we would have $g^{j+1}=g^{j}$ implying, in view of Eq. (6.21), that $g^{j}=0$. Therefore, we have using Eqs. (6.21) and (6.22)

$$
g^{i^{\prime}} Q d^{j}=\frac{1}{\alpha^{j}} g^{i^{\prime}}\left(g^{j+1}-g^{j}\right)=\left\{\begin{array}{lll}
0 & \text { if } & j=0, \ldots, i-2 \\
\frac{1}{\alpha^{j}} g^{i^{\prime}} g^{i} & \text { if } & j=i-1
\end{array}\right.
$$

and also that

$$
d^{j^{\prime}} Q d^{j}=\frac{1}{\alpha^{j}} d^{j^{\prime}}\left(g^{j+1}-g^{j}\right)
$$

Substituting the last two relations in the Gram-Schmidt formula (6.18) we obtain

$$
d^{k}=-g^{k}+\beta^{k} d^{k-1}
$$

where

$$
\beta^{k}=\frac{g^{k^{\prime}} g^{k}}{d^{k-1^{\prime}}\left(g^{k}-g^{k-1}\right)}
$$

From Eq. (6.23) we have $d^{k-1}=-g^{k-1}+\beta^{k-1} d^{k-2}$. Using this equation, and the orthogonality of $g^{k}$ and $g^{k-1}$, and of $d^{k-2}$ and $g^{k}-g^{k-1}$ [cf. Eqs. (6.20) and (6.21)], the denominator in Eq. (6.24) is written as $g^{k-1^{\prime}} g^{k-1}$, and the desired formula for $\beta^{k}$ follows. Q.E.D.

Note that by using the orthogonality of $g^{k}$ and $g^{k-1}$ the formula

$$
\beta^{k}=\frac{g^{k^{\prime}} g^{k}}{g^{k-1^{\prime}} g^{k-1}}
$$

---

of Prop. 1.6.1 can also be written as

$$
\beta^{k}=\frac{g^{k^{\prime}}\left(g^{k}-g^{k-1}\right)}{g^{k-1^{\prime}} g^{k-1}}
$$

While the alternative formulas (6.25), and (6.26) produce the same results for quadratic problems, their differences become significant when the conjugate gradient method is extended to nonquadratic problems, as we will discuss shortly.

# Preconditioned Conjugate Gradient Method 

This method is really the conjugate gradient method implemented in a new coordinate system. Suppose we make a change of variables, $x=S y$, where $S$ is an invertible symmetric $n \times n$ matrix, and we apply the conjugate gradient method to the equivalent problem

$$
\begin{array}{ll}
\operatorname{minimize} & h(y)=f(S y)=\frac{1}{2} y^{\prime} S Q S y-b^{\prime} S y \\
\text { subject to } & y \in \Re^{n}
\end{array}
$$

The method is described by

$$
y^{k+1}=y^{k}+\alpha^{k} \tilde{d}^{k}
$$

where $\alpha^{k}$ is obtained by line minimization and $\tilde{d}^{k}$ is generated by [cf. Eqs. (6.23) and (6.25)]

$$
\tilde{d}^{0}=-\nabla h\left(y^{0}\right), \quad \tilde{d}^{k}=-\nabla h\left(y^{k}\right)+\beta^{k} \tilde{d}^{k-1}, \quad k=1, \ldots, n-1
$$

where

$$
\beta^{k}=\frac{\nabla h\left(y^{k}\right)^{\prime} \nabla h\left(y^{k}\right)}{\nabla h\left(y^{k-1}\right)^{\prime} \nabla h\left(y^{k-1}\right)}
$$

Setting $x^{k}=S y^{k}, \nabla h\left(y^{k}\right)=S g^{k}, d^{k}=S \tilde{d}^{k}$, and $H=S^{2}$, we obtain from Eqs. (6.27)-(6.29) the equivalent method

$$
\begin{gathered}
x^{k+1}=x^{k}+\alpha^{k} d^{k} \\
d^{0}=-H g^{0}, \quad d^{k}=-H g^{k}+\beta^{k} d^{k-1}, \quad k=1, \ldots, n-1
\end{gathered}
$$

where

$$
\beta^{k}=\frac{g^{k^{\prime}} H g^{k}}{g^{k-1^{\prime}} H g^{k-1}}
$$

and $\alpha^{k}$ is obtained by line minimization.
The method described by the above equations is called the preconditioned conjugate gradient method with scaling matrix $H$. To see that this

---

method is a conjugate direction method, note that since $\nabla^{2} h(y)=S Q S$, the vectors $d^{0} \ldots d^{n-1}$ are ( $S Q S$ )-conjugate. Since $d^{k}=S d^{k}$, we obtain that $d^{0}, \ldots, d^{n-1}$ are $Q$-conjugate. Therefore, the scaled method terminates with the minimum of $f$ after at most $n$ iterations, just as the ordinary conjugate gradient method. The motivation for scaling is to improve the rate of convergence within an $n$-iteration cycle (see the following analysis). This is important for a nonquadratic problem, but it may be important even for a quadratic problem if $n$ is large and we want to obtain an approximate solution without waiting for the method to terminate.

# Application to Nonquadratic Problems 

The conjugate gradient method can be applied to the nonquadratic problem

$$
\begin{array}{ll}
\operatorname{minimize} & f(x) \\
\text { subject to } & x \in \Re^{n}
\end{array}
$$

It takes the form

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}
$$

where $\alpha^{k}$ is obtained by line minimization

$$
f\left(x^{k}+\alpha^{k} d^{k}\right)=\min _{\alpha} f\left(x^{k}+\alpha d^{k}\right)
$$

and $d^{k}$ is generated by

$$
d^{k}=-\nabla f\left(x^{k}\right)+\beta^{k} d^{k-1}
$$

The most common way to compute $\beta^{k}$ is [cf. Eq. (6.26)]

$$
d^{k}=\frac{\nabla f\left(x^{k}\right)^{\prime}\left(\nabla f\left(x^{k}\right)-\nabla f\left(x^{k-1}\right)\right)}{\nabla f\left(x^{k-1}\right)^{\prime} \nabla f\left(x^{k-1}\right)}
$$

The direction $d^{k}$ generated by the formula $d^{k}=-\nabla f\left(x^{k}\right)+\beta^{k} d^{k-1}$ is a direction of descent, since from Eq. (6.31) we obtain $\nabla f\left(x^{k}\right)^{\prime} d^{k-1}=0$, so that

$$
\nabla f\left(x^{k}\right)^{\prime} d^{k}=-\left\|\nabla f\left(x^{k}\right)\right\|^{2}+\beta^{k} \nabla f\left(x^{k}\right)^{\prime} d^{k-1}=-\left\|\nabla f\left(x^{k}\right)\right\|^{2}
$$

For nonquadratic problems, the formula (6.36) is typically superior to alternative formulas such as

$$
d^{k}=\frac{\nabla f\left(x^{k}\right)^{\prime} \nabla f\left(x^{k}\right)}{\nabla f\left(x^{k-1}\right)^{\prime} \nabla f\left(x^{k-1}\right)}
$$

2. 6.25 A hernstic explanation is that due to nonquadratic terms in the obiective function and possibly inaccurate line searches, conjugacy

---

of the generated directions is progressively lost and a situation may arise where the method "jams" in the sense that the generated direction $d^{k}$ is nearly orthogonal to the gradient $\nabla f\left(x^{k}\right)$. When this occurs, we have $\nabla f\left(x^{k+1}\right) \simeq \nabla f\left(x^{k}\right)$. In that case, the scalar $\beta^{k+1}$, generated by

$$
\beta^{k+1}=\frac{\nabla f\left(x^{k+1}\right)^{\prime}\left(\nabla f\left(x^{k+1}\right)-\nabla f\left(x^{k}\right)\right)}{\nabla f\left(x^{k}\right)^{\prime} \nabla f\left(x^{k}\right)}
$$

will be nearly zero and the next direction $d^{k+1}=-\nabla f\left(x^{k+1}\right)+\beta^{k+1} d^{k}$ will be close to $-\nabla f\left(x^{k+1}\right)$ thereby breaking the jam. By contrast, when Eq. (6.37) is used, under the same circumstances the method typically continues to jam.

Regardless of the direction update formula used, one must deal with the loss of conjugacy that results from nonquadratic terms in the cost function. The conjugate gradient method is often employed in problems where the number of variables $n$ is large, and it is not unusual for the method to start generating nonsensical and inefficient directions of search after a few iterations. For this reason it is important to operate the method in cycles of conjugate direction steps, with the first step in the cycle being a steepest descent step. Some possible restarting policies are:
(a) Restart with a steepest descent step $n$ iterations after the preceding restart.
(b) Restart with a steepest descent step $k$ iterations after the preceding restart with $k<n$. This is recommended when the problem has special structure so that the resulting method has good convergence rate (see the following Prop. 1.6.2).
(c) Restart with a steepest descent step if either $n$ iterations have taken place since the preceding restart or if

$$
\left|\nabla f\left(x^{k}\right)^{\prime} \nabla f\left(x^{k-1}\right)\right|>\gamma\left\|\nabla f\left(x^{k-1}\right)\right\|^{2}
$$

where $\gamma$ is a fixed scalar with $0<\gamma<1$. The above relation is a test on loss of conjugacy, for if the generated directions were conjugate then we would have $\nabla f\left(x^{k}\right)^{\prime} \nabla f\left(x^{k-1}\right)=0$.

Note that in all these restart procedures the steepest descent iteration serves as a spacer step and guarantees global convergence (Prop. 1.2.6 in Section 1.2). If the scaled version of the conjugate gradient method is used, then a scaled steepest descent iteration is used to restart a cycle. The scaling matrix may change at the beginning of a cycle but should remain unchanged during the cycle.

An important practical issue relates to the line search accuracy that is necessary for efficient computation. On one hand, an accurate line search is needed to limit the loss of direction conjugacy and the attendant deterioration of convergence rate. On the other hand, insisting on a very accurate

---

line search can be computationally expensive. Some trial and error may therefore be required in practice. For a discussion of implementations that are tolerant of line search inaccuracies see [Per78] and [Sha78]. For a computational study comparing different implementations, see [PaG86].

# Conjugate Gradient-Like Methods for Linear Systems* 

The conjugate gradient method can be used to solve the linear system of equations

$$
A x=b
$$

where $A$ is an invertible $n \times n$ matrix and $b$ is a given vector in $\Re^{n}$. One way to do this is to apply the conjugate gradient method to the positive definite quadratic optimization problem

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2} x^{\prime} A^{\prime} A x-b^{\prime} A x \\
& \text { subject to } x \in \Re^{n}
\end{aligned}
$$

which corresponds to the equivalent linear system $A^{\prime} A x=A^{\prime} b$. This, however, has several disadvantages, including the need to form the matrix $A^{\prime} A$, which may have a much less favorable sparsity structure than $A$.

An alternative possibility is to introduce the vector $z$ defined by

$$
x=A^{\prime} z
$$

and to solve the system $A A^{\prime} z=b$ or equivalently the positive definite quadratic problem

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2} z^{\prime} A A^{\prime} z-b^{\prime} z \\
& \text { subject to } z \in \Re^{n}
\end{aligned}
$$

whose cost function gradient is zero at $z$ if and only if $A A^{\prime} z=b$. By streamlining the computations, it is possible to write the conjugate gradient method for the preceding problem directly in terms of the vector $x$, and without explicitly forming the product $A A^{\prime}$. The resulting method is known as Craig's method; it is given by the following iteration where $H$ is a positive definite symmetric preconditioning matrix:

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}, \quad \alpha^{k}=\frac{r^{k^{\prime}} r^{k}}{d^{k^{\prime}} d^{k}}
$$

where the vectors $r^{k}$ and $d^{k}$ are generated by the recursions

$$
r^{k+1}=r^{k}+\alpha^{k} H A d^{k}, \quad d^{k+1}=-A^{\prime} H r^{k+1}+\frac{r^{k+1^{\prime}} r^{k+1}}{r^{k^{\prime}} r^{k}} d^{k}
$$

with the initial conditions

$$
r^{0}=H\left(A x^{0}-b\right), \quad d^{0}=-A^{\prime} H r^{0}
$$

---

The verification of these equations is left for the reader.
There are other conjugate gradient-like methods for the system $A x=$ $b$, which are not really equivalent to the conjugate gradient method for any quadratic optimization problem. One possibility, due to [SaS86], known as the Generalized Minimum Residual method (GMRES), is to start with a vector $x^{0}$ and obtain $x^{k}$ as the vector that minimizes $\|A x-b\|^{2}$ over the linear manifold $x^{0}+S^{k}$, where

$$
S^{k}=\left(\text { subspace spanned by the vectors } r, A r, A^{2} r, \ldots, A^{k-1} r\right)
$$

and $r$ is the initial residual

$$
r=A x^{0}-b
$$

This successive subspace minimization process can be efficiently implemented, but we will not get into the details further (see [SaS86]). It can be shown that $x^{k}$ is a solution of the system $A x=b$ if and only if $A^{k} r$ belongs to the subspace $S^{k}$ (write the minimization of $\|A x-b\|^{2}$ over $x^{0}+S^{k}$ as the equivalent minimization of $|\xi-r|^{2}$ over all $\xi$ in the subspace $A S^{k}$ ). Thus if none of the vectors $x^{0}, \ldots, x^{n-2}$ is a solution, the subspace $S^{n-1}$ is equal to $\Re^{n}$, implying that $x^{n-1}$ is an unconstrained minimum of $\|A x-b\|^{2}$, and therefore solves the system $A x=b$. It follows that the method will terminate after at most $n$ iterations.

GMRES can be viewed as a conjugate gradient method only in the special case where $A$ is positive definite and symmetric. In that case it can be shown that the method is equivalent to a preconditioned conjugate gradient method applied to the quadratic cost $\|A x-b\|^{2}$. This is based on the expanding subspace minimization property of the conjugate gradient method (see also Exercise 6.4). Note, however, that GMRES can be used for any matrix $A$ that is invertible.

# Rate of Convergence of the Conjugate Gradient Method* 

There are a number of convergence rate results for the conjugate gradient method. Since the method terminates in at most $n$ steps for a quadratic cost, one would expect that when viewed in cycles of $n$ steps, its rate of convergence for a nonquadratic cost would be comparable to the rate of Newton's method. Indeed there are results which roughly state that if the method is restarted every $n$ iterations and $\left\{x^{k}\right\}$ converges to a nonsingular local minimum $x^{*}$, then the error $e^{k}=\left\|x^{n k}-x^{*}\right\|$ converges superlinearly. (Note that here the error is considered at the end of cycles of $n$ iterations rather than at the end of each iteration.) Such results are reassuring but not terribly interesting because the conjugate gradient method is most useful in problems where $n$ is large (see the discussion at the end of Section 1.7), and for such problems, one hopes that practical

---

convergence will occur after fewer than $n$ iterations. Therefore, the singlestep rate of convergence of the method is more interesting than its rate of convergence in terms of $n$-step cycles. The following analysis gives a result of this type, based on an interpretation of the conjugate gradient method as an optimal process.

Assume that the cost is positive definite quadratic of the form

$$
f(x)=\frac{1}{2} x^{\prime} Q x
$$

(To simplify the following exposition, we have assumed that the linear term $b^{\prime} x$ is zero, but with minor modifications, the following analysis holds also when $b \neq 0$.) Let $g^{i}$ denote as usual the gradient $\nabla f\left(x^{i}\right)$ and consider an algorithm of the form

$$
\begin{gathered}
x^{1}=x^{0}+\gamma^{00} g^{0} \\
x^{2}=x^{0}+\gamma^{10} g^{0}+\gamma^{11} g^{1} \\
\ldots \\
x^{k+1}=x^{0}+\gamma^{k 0} g^{0}+\cdots+\gamma^{k k} g^{k}
\end{gathered}
$$

where $\gamma^{i j}$ are arbitrary scalars. Since $g^{i}=Q x^{i}$, we see that for suitable scalars $c^{k i}$, the above algorithm can be written for all $k$ as

$$
x^{k+1}=x^{0}+c^{k 0} Q x^{0}+c^{k 1} Q^{2} x^{0}+\cdots+c^{k k} Q^{k+1} x^{0}=\left(I+Q P^{k}(Q)\right) x^{0}
$$

where $P^{k}$ is a polynomial of degree $k$.
Among algorithms of the form (6.39), the conjugate gradient method is optimal in the sense that for every $k$, it minimizes $f\left(x^{k+1}\right)$ over all sets of coefficients $\gamma^{k 0}, \ldots, \gamma^{k k}$. It follows from the equation above that in the conjugate gradient method we have, for every $k$,

$$
f\left(x^{k+1}\right)=\min _{P^{k}} \frac{1}{2} x^{0^{\prime}} Q\left(I+Q P^{k}(Q)\right)^{2} x^{0}
$$

Let $\lambda_{1}, \ldots, \lambda_{n}$ be the eigenvalues of $Q$, and let $e_{1}, \ldots, e_{n}$ be corresponding orthogonal eigenvectors, normalized so that $\left\|e_{i}\right\|=1$. Since $e_{1}, \ldots, e_{n}$ form a basis, any vector $x^{0} \in \Re^{n}$ can be written as

$$
x^{0}=\sum_{i=1}^{n} \xi_{i} e_{i}
$$

for some scalars $\xi_{i}$. Since

$$
Q x^{0}=\sum_{i=1}^{n} \xi_{i} Q e_{i}=\sum_{i=1}^{n} \xi_{i} \lambda_{i} e_{i}
$$

---

we have, using the orthogonality of $e_{1}, \ldots, e_{n}$ and the fact $\left\|e_{i}\right\|=1$,

$$
f\left(x^{0}\right)=\frac{1}{2} x^{0^{\prime}} Q x^{0}=\frac{1}{2}\left(\sum_{i=1}^{n} \xi_{i} c_{i}\right)^{\prime}\left(\sum_{i=1}^{n} \xi_{i} \lambda_{i} c_{i}\right)=\frac{1}{2} \sum_{i=1}^{n} \lambda_{i} \xi_{i}^{2}
$$

Applying the same process to E(. (6.40), we obtain for any polynomial $P^{k}$ of degree $k$

$$
f\left(x^{k+1}\right) \leq \frac{1}{2} \sum_{i=1}^{n}\left(1+\lambda_{i} P^{k}\left(\lambda_{i}\right)\right)^{2} \lambda_{i} \xi_{i}^{2}
$$

and it follows that

$$
f\left(x^{k+1}\right) \leq \max _{i}\left(1+\lambda_{i} P^{k}\left(\lambda_{i}\right)\right)^{2} f\left(x^{0}\right), \quad \forall P^{k}, k
$$

One can use this relationship for different choices of polynomials $P^{k}$ to obtain a number of convergence rate results. We provide one such result, which shows that the first $k$ conjugate gradient iterations in an $n$-iteration cycle eliminate the effect of the $k$ largest eigenvalues of $Q$.

Proposition 1.6.2: Assume that $Q$ has $n-k$ eigenvalues in an interval $[a, b]$ with $a>0$, and the remaining $k$ eigenvalues are greater than $b$. Then for every $x^{0}$, the vector $x^{k+1}$ generated after $k+1$ steps of the conjugate gradient method satisfies

$$
f\left(x^{k+1}\right) \leq\left(\frac{b-a}{b+a}\right)^{2} f\left(x^{0}\right)
$$

This relation also holds for the preconditioned conjugate gradient method (6.30)-(6.32) if the eigenvalues of $Q$ are replaced by those of $H^{1 / 2} Q H^{1 / 2}$.

Proof: Let $\lambda_{1}, \ldots \lambda_{k}$ be the eigenvalues of $Q$ that are greater than $b$ and consider the polynomial $P^{k}$ defined by

$$
1+\lambda P^{k}(\lambda)=\frac{2}{(a+b) \lambda_{1} \cdots \lambda_{k}}\left(\frac{a+b}{2}-\lambda\right)\left(\lambda_{1}-\lambda\right) \cdots\left(\lambda_{k}-\lambda\right)
$$

Since $1+\lambda_{i} P^{k}\left(\lambda_{i}\right)=0$, we have, using Eqs. (6.41), (6.43), and a simple calculation,

$$
f\left(x^{k+1}\right) \leq \max _{a \leq \lambda \leq b} \frac{\left(\lambda-\frac{1}{2}(a+b)\right)^{2}}{\left(\frac{1}{2}(a+b)\right)^{2}} f\left(x^{0}\right)=\left(\frac{b-a}{b+a}\right)^{2} f\left(x^{0}\right)
$$

Q.E.D.

---

One consequence of the above proposition is that if the eigenvalues of $Q$ take only $k$ distinct values then the conjugate gradient method will find the minimum of the quadratic function $f$ in at most $k$ iterations. (Take $a=b$.) Some other possibilities are explored in the problem section.

It is worth mentioning two more rate of convergence results regarding the conjugate gradient method as applied to the positive definite quadratic function

$$
f(x)=\frac{1}{2}\left(x-x^{*}\right)^{\prime} Q\left(x-x^{*}\right)
$$

Let $M$ and $m$ be the largest and smallest eigenvalues of $f$, respectively. Then, for any starting point $x^{0}$ and any iteration index $k$, it can be shown (see [Pol87]) that

$$
\begin{gathered}
\left\|x^{k}-x^{*}\right\| \leq 2\left(\frac{M}{m}\right)^{1 / 2}\left(\frac{\sqrt{M}-\sqrt{m}}{\sqrt{M}+\sqrt{m}}\right)^{k}\left\|x^{0}-x^{*}\right\| \\
f\left(x^{k}\right) \leq \frac{M\left\|x^{0}-x^{*}\right\|^{2}}{2(2 k+1)^{2}}
\end{gathered}
$$

These relations again suggest a more favorable convergence rate than the one of steepest descent; compare with the results of Section 1.3.

# E XERCISES 

## 6.1

Show that the Gram-Schmidt procedure has the projection property stated in Fig. 1.6.2.

## 6.2 [Ber74]

Let $Q$ have the form

$$
Q=M+\sum_{i=1}^{k} v_{i} v_{i}^{\prime}
$$

where $M$ is positive definite, and $v_{i}$ are some vectors in $\Re^{n}$. Show that the vector $x^{k+1}$ generated after $k+1$ steps of the conjugate gradient method satisfies

$$
f\left(x^{k+1}\right) \leq\left(\frac{b-a}{b+a}\right)^{2} f\left(x^{0}\right)
$$

where $a$ and $b$ are the smallest and largest eigenvalues of $M$, respectively. Show also that the vector $x^{k+1}$ generated by the preconditioned conjugate gradient method with $H=M^{-1}$ minimizes $f$. Hint: Use the interlocking eigenvalues lemma [Prop. A.18(d) in Appendix A].

---

# 6.3 (Hessian with Clustered Eigenvalues [Ber82a]) 

Assume that $Q$ has all its eigenvalues concentrated at $k$ intervals of the form

$$
\left[z_{i}-\delta_{i}, z_{i}+\delta_{i}\right], \quad i=1, \ldots, k
$$

where we assume that $\delta_{i} \geq 0, i=1, \ldots, k, 0<z_{1}-\delta_{1}$, and

$$
0<z_{1}<z_{2}<\cdots<z_{k}, \quad z_{i}+\delta_{i} \leq z_{i+1}-\delta_{i+1}, \quad i=1, \ldots, k-1
$$

Show that the vector $x^{k-1}$ generated after $k+1$ steps of the conjugate gradient method satisfies

$$
f\left(x^{k+1}\right) \leq R f\left(x^{0}\right)
$$

where

$$
\begin{aligned}
R=\max & \left\{\frac{\delta_{k}^{2}}{z_{1}^{2}}: \frac{\delta_{2}^{2}\left(z_{2}+\delta_{2}-z_{1}\right)^{2}}{z_{1}^{2} z_{2}^{2}}\right. \\
& \left.\frac{\delta_{1}^{2}\left(z_{k}+\delta_{k}-z_{1}\right)^{2}\left(z_{k}+\delta_{k}-z_{2}\right)^{2} \cdots\left(z_{k}+\delta_{k}-z_{k-1}\right)^{2}}{z_{1}^{2} z_{2}^{2} \cdots z_{k}^{2}}\right\}
\end{aligned}
$$

6.4

Consider the conjugate gradient method applied to the minimization of $f(x)=$ $\frac{1}{2} x^{2} Q x-h^{2} x$, where $Q$ is positive definite and symmetric. Show that the iterate $x^{k}$ minimizes $f$ over the linear manifold

$$
x^{0}+\left(\text { subspace spanned by } g^{0}, Q g^{0}, \ldots, Q^{k-1} g^{0}\right)
$$

where $g^{0}=\nabla f\left(x^{0}\right)$
6.5

Let $f$ be positive definite quadratic. Consider the following method: The first iteration is a steepest descent iteration with the stepsize determined by line minimization. For $k=2, \ldots, n$, the $k$ th iteration finds $x^{k}$ that minimizes $f$ over the two-dimensional subspace spanned by $g^{k}$ and $x^{k}-x^{k-1}$. Show that this method is equivalent to the conjugate gradient method.
6.6

Suppose that $d^{0}, \ldots, d^{k}$ are $Q$-conjugate directions, let $x^{1}, \ldots, x^{k+1}$ be the vectors generated by the corresponding conjugate direction method, and assume that $x^{i+1} \neq x^{i}$ for all $i=0, \ldots, k$. Show that a vector $d^{k+1}$ is $Q$ conjugate to $d^{0}, \ldots, d^{k}$ if and only if $d^{k-1} \neq 0$ and $d^{k-1}$ is orthogonal to the gradient differences $g^{i+1}-g^{i}, i=0, \ldots, k$.

---

6.7

Describe the behavior of the conjugate gradient method for a positive semidefinite quadratic function. Consider the case where there is no optimal solution and the case where there are infinitely many optimal solutions.
6.8

Let $f(x)=\frac{1}{2} x^{\prime} Q x-b^{\prime} x$, where $Q$ is positive definite and symmetric. Suppose that $x_{1}$ and $x_{2}$ minimize $f$ over linear manifolds that are parallel to subspaces $S_{1}$ and $S_{2}$, respectively. Show that if $x_{1} \neq x_{2}$, then $x_{1}-x_{2}$ is $Q$-conjugate to all vectors in the intersection of $S_{1}$ and $S_{2}$. Use this property to construct a conjugate direction method that does not evaluate gradients and uses only line minimizations.

# 1.7 QUASI-NEWTON METHODS 

Quasi-Newton methods are gradient methods of the form

$$
\begin{aligned}
& x^{k+1}=x^{k}+\alpha^{k} d^{k} \\
& d^{k}=-D^{k} \nabla f\left(x^{k}\right)
\end{aligned}
$$

where $D^{k}$ is a positive definite matrix, which may be adjusted from one iteration to the next so that the direction $d^{k}$ tends to approximate the Newton direction. Some of these methods are quite popular because they typically converge fast, while avoiding the second derivative calculations associated with Newton's method. Their main drawback relative to the conjugate gradient method is that they require storage of the matrix $D^{k}$ as well as the matrix-vector multiplication overhead associated with the calculation of the direction $d^{k}$ (see the subsequent discussion).

An important idea for many quasi-Newton methods is that two successive iterates $x^{k}, x^{k+1}$ together with the corresponding gradients $\nabla f\left(x^{k}\right)$, $\nabla f\left(x^{k+1}\right)$, yield curvature information by means of the approximate relation

$$
q^{k} \approx \nabla^{2} f\left(x^{k+1}\right) p^{k}
$$

where

$$
\begin{gathered}
p^{k}=x^{k+1}-x^{k} \\
q^{k}=\nabla f\left(x^{k+1}\right)-\nabla f\left(x^{k}\right)
\end{gathered}
$$

---

In particular, given $n$ linearly independent iteration increments $p^{0}, \ldots, p^{n-1}$ together with the corresponding gradient increments $q^{0}, \ldots, q^{n-1}$, we can obtain approximately the Hessian as

$$
\nabla^{2} f\left(x^{n}\right) \approx\left[q^{0} \cdots q^{n-1}\right]\left[p^{0} \cdots p^{n-1}\right]^{-1}
$$

and the inverse Hessian as

$$
\nabla^{2} f\left(x^{n}\right)^{-1} \approx\left[p^{0} \cdots p^{n-1}\right]\left[q^{0} \cdots q^{n-1}\right]^{-1}
$$

When the cost is quadratic, this relation is exact. Many interesting quasiNewton methods use similar but more sophisticated ways to build curvature information into the matrix $D^{k}$ so that it progressively approaches the inverse Hessian.

In the most popular class of quasi-Newton methods, the matrix $D^{k+1}$ is obtained from $D^{k}$, and the vectors $p^{k}$ and $q^{k}$ by means of the equation

$$
D^{k+1}=D^{k}+\frac{p^{k} p^{k^{\prime}}}{p^{k^{\prime}} q^{k}}-\frac{D^{k} q^{k} q^{k^{\prime}} D^{k}}{q^{k^{\prime}} D^{k} q^{k}}+\xi^{k} \tau^{k} v^{k} v^{k^{\prime}}
$$

where

$$
\begin{gathered}
v^{k}=\frac{p^{k}}{p^{k^{\prime}} q^{k}}-\frac{D^{k} q^{k}}{\tau^{k}} \\
\tau^{k}=q^{k^{\prime}} D^{k} q^{k}
\end{gathered}
$$

the scalars $\xi^{k}$ satisfy, for all $k$,

$$
0 \leq \xi^{k} \leq 1
$$

and $D^{0}$ is an arbitrary positive definite matrix. The scalars $\xi^{k}$ parameterize the method. If $\xi^{k}=0$ for all $k$, we obtain the Davidon-Fletcher-Powell (DFP) method, which is historically the first quasi-Newton method. If $\xi^{k}=1$ for all $k$, we obtain the Broyden-Fletcher-Goldfarb-Shanno (BFGS) method, for which there is substantial evidence that it is the best general purpose quasi-Newton method currently known.

We first show that under a mild assumption, the matrices $D^{k}$ generated by Eq. (7.6) are positive definite. This is a very important property, since it guarantees that $d^{k}$ is a descent direction.

Proposition 1.7.1: If $D^{k}$ is positive definite and the stepsize $\alpha^{k}$ is chosen so that $x^{k+1}$ satisfies

$$
\nabla f\left(x^{k}\right)^{\prime} d^{k}<\nabla f\left(x^{k+1}\right)^{\prime} d^{k}
$$

then $D^{k+1}$ as given by Eq. (7.6) is positive definite.

---

Note: If $x^{k}$ is not a stationary point, we have $\nabla f\left(x^{k}\right)^{\prime} d^{k}<0$, so in order to satisfy condition (7.10), it is sufficient to carry out the line search to a point where

$$
\left|\nabla f\left(x^{k+1}\right)^{\prime} d^{k}\right|<\left|\nabla f\left(x^{k}\right)^{\prime} d^{k}\right|
$$

In particular, if $a^{k}$ is determined by the line minimization rule, then we have $\nabla f\left(x^{k+1}\right)^{\prime} d^{k}=0$ and Eq. (7.10) is satisfied.

Proof: We first note that Eq. (7.10) implies that $\alpha^{k} \neq 0, q^{k} \neq 0$, and

$$
p^{k^{\prime}} q^{k}=\alpha^{k} d^{k^{\prime}}\left(\nabla f\left(x^{k+1}\right)-\nabla f\left(x^{k}\right)\right)>0
$$

Thus all denominator terms in Eqs. (7.6) and (7.7) are nonzero, and $D^{k+1}$ is well defined.

Now for any $z \neq 0$ we have

$$
z^{\prime} D^{k+1} z=z^{\prime} D^{k} z+\frac{\left(z^{\prime} p^{k}\right)^{2}}{p^{k^{\prime}} q^{k}}-\frac{\left(q^{k^{\prime}} D^{k} z\right)^{2}}{q^{k^{\prime}} D^{k} q^{k}}+\xi^{k} \tau^{k}\left(v^{k^{\prime}} z\right)^{2}
$$

Using the notation $a=\left(D^{k}\right)^{12} z, b=\left(D^{k}\right)^{12} q^{k}$, this equation is written as

$$
z^{\prime} D^{k+1} z=\frac{\|a\|^{2}\|b\|^{2}-\left(a^{\prime} b\right)^{2}}{\|b\|^{2}}+\frac{\left(z^{\prime} p^{k}\right)^{2}}{p^{k^{\prime}} q^{k}}+\xi^{k} \tau^{k}\left(v^{k^{\prime}} z\right)^{2}
$$

From Eqs. (7.8), and (7.11), and the Schwartz inequality [Eq. (A.2) in Appendix A], we have that all terms on the right-hand side of Eq. (7.13) are nonnegative. In order that $z^{\prime} D^{k+1} z>0$, it will suffice to show that we cannot have simultaneously

$$
\|a\|^{2}\|b\|^{2}=\left(a^{\prime} b\right)^{2} \quad \text { and } \quad z^{\prime} p^{k}=0
$$

Indeed if $\left\|a\right\|^{2}\|b\|^{2}=\left(a^{\prime} b\right)^{2}$, we must have $a=\lambda b$ or equivalently, $z=\lambda q^{k}$. Since $z \neq 0$, it follows that $\lambda \neq 0$, so if $z^{\prime} p^{k}=0$, we must have $q^{k^{\prime}} p^{k}=0$, which is impossible by Eq. (7.11). Q.E.D.

An important property of the algorithm is that when applied to the positive definite quadratic function $f(x)=\frac{1}{2} x^{\prime} Q x-b^{\prime} x$, with the stepsize $a^{k}$ determined by line minimization, it generates a $Q$-conjugate direction sequence, while simultaneously constructing the inverse Hessian $Q^{-1}$ after $n$ iterations. This is the subject of the next proposition.

---

Proposition 1.7.2: Let $\left\{x^{k}\right\},\left\{d^{k}\right\}$, and $\left\{D^{k}\right\}$ be sequences jponeu ated by the Quasi-Newton algorithm (7.1)-(7.2). (7.6)-(7.9). applied to minimization of the positive definite quadratic function

$$
f(x)=\frac{1}{2} x^{\prime} Q x-b^{\prime} x
$$

with $\alpha^{k}$ chosen by

$$
f\left(x^{k}+\alpha^{k} d^{k}\right)=\min _{i} f\left(x^{k}+\alpha d^{k}\right)
$$

Assume that none of the vectors $x^{0}, \ldots, x^{n-1}$ is optimal. Then:
(a) The vectors $d^{0}, \ldots, d^{n-1}$ are $Q$-conjugate.
(b) There holds

$$
D^{n}=Q^{-1}
$$

Proof: We will show that for all $k$

$$
\begin{gathered}
d^{i^{\prime}} Q d^{j}=0, \quad 0 \leq i<j \leq k \\
D^{k+1} Q p^{i}=p^{i}, \quad 0 \leq i \leq k
\end{gathered}
$$

Equation (7.15) proves part (a) and it can be shown that Eq. (7.16) proves part (b). Indeed, since for $i<n$, none of the vectors $x^{i}$ is optimal and $d^{i}$ is a descent direction (cf. Eq. (7.2) and Prop. 1.7.1), we have that $p^{i} \neq 0$. Since $p^{i}=\alpha^{i} d^{i}$ and $d^{0}, \ldots, d^{n-1}$ are $Q$-conjugate, it follows that $p^{0}, \ldots, p^{n-1}$ are linearly independent and therefore. Eq. (7.16) implies that $D^{n} Q$ is equal to the identity matrix.

We first verify that

$$
D^{k+1} Q p^{k}=p^{k}, \quad \forall k
$$

From the equation $Q p^{k}=q^{k}$ and the updating formula (7.6), we have

$$
\begin{aligned}
D^{k+1} Q p^{k} & =D^{k+1} q^{k}=D^{k} q^{k}+\frac{p^{k} p^{k^{\prime}} q^{k}}{p^{k^{\prime}} q^{k}}-\frac{D^{k} q^{k} q^{k^{\prime}} D^{k} q^{k}}{q^{k^{\prime}} D^{k} q^{k}}+\xi^{k} \tau^{k} v^{k} v^{k^{\prime}} q^{k} \\
& =p^{k}+\xi^{k} \tau^{k} v^{k} v_{k^{\prime}} q^{k}
\end{aligned}
$$

From Eqs. (7.7) and (7.8), we have that $r^{k^{\prime}} q^{k}=0$ and Eq. (7.17) follows.
We now show Eqs. (7.15) and (7.16) simultaneously by induction. For $k=0$ there is nothing to show for Eq. (7.15), while Eq. (7.16) holds in view of Eq. (7.17). Assuming that Eqs. (7.15) and (7.16) hold for $k$, we prove them for $k+1$. We have, for $i<k$,

$$
\nabla f\left(x^{k+1}\right)=\nabla f\left(x^{i+1}\right)+Q\left(p^{i+1}+\cdots+p^{k}\right)
$$

---

The vector $p^{i}$ is orthogonal to each vector in the right-hand side of this equation; it is orthogonal to $Q p^{i+1}, \ldots, Q p^{k}$ because $p^{0}, \ldots, p^{k}$ are $Q$ conjugate (since $p^{i}=\alpha^{i} d^{i}$ ) and it is orthogonal to $\nabla f\left(x^{i+1}\right)$ because of the line minimization property of the stepsize [cf. Eq. (7.14)]. Therefore from Eq. (7.18) we obtain

$$
p^{i^{\prime}} \nabla f\left(x^{k+1}\right)=0, \quad 0 \leq i<k
$$

From this equation and Eq. (7.16),

$$
p^{i^{\prime}} Q D^{k+1} \nabla f\left(x^{k+1}\right)=0, \quad 0 \leq i \leq k
$$

and since $p^{i}=\alpha^{i} d^{i}$ and $d^{k+1}=-D^{k+1} \nabla f\left(x^{k+1}\right)$, we obtain

$$
d^{i^{\prime}} Q d^{k+1}=0, \quad 0 \leq i \leq k
$$

This proves Eq. (7.15) for $k+1$.
From the induction hypothesis (7.16) and Eq. (7.19), we have for all $i$ with $0 \leq i \leq k$,

$$
q^{k+1^{\prime}} D^{k+1} Q p^{i}=q^{k+1^{\prime}} p^{i}=p^{k+1^{\prime}} Q p^{i}=\alpha^{k+1} \alpha^{i} d^{k+1^{\prime}} Q d^{i}=0
$$

From Eq. (7.6), we have, for $0 \leq i \leq k$,

$$
\begin{aligned}
D^{k+2} q^{i}= & D^{k+1} q^{i}+\frac{p^{k+1} p^{k+1^{\prime}} q^{i}}{p^{k+1^{\prime}} q^{k+1}}-\frac{D^{k+1} q^{k+1} q^{k+1^{\prime}} D^{k+1} q^{i}}{q^{k+1^{\prime}} D^{k+1} q^{k+1}} \\
& +\xi^{k+1} \tau^{k+1} v^{k+1} v^{k+1^{\prime}} q^{i}
\end{aligned}
$$

Since $p^{k+1^{\prime}} q^{i}=p^{k+1^{\prime}} Q p^{i}=\alpha^{k+1} \alpha^{i} d^{k+1^{\prime}} Q d^{i}=0$, we see that the second term in the right-hand side of Eq. (7.21) is zero. Similarly, Eq. (7.16) implies that $q^{k+1^{\prime}} D^{k+1} q^{i}=q^{k+1^{\prime}} D^{k+1} Q p^{i}=q^{k+1^{\prime}} p^{i}=p^{k+1^{\prime}} Q p^{i}=0$ and we see that the third term in the right-hand side of Eq. (7.21) is zero. Finally, a similar argument using the definition (7.7) of $v^{k+1}$, shows that the fourth term in the right-hand side of Eq. (7.21) is zero as well. Therefore, Eqs. (7.21) and (7.16) yield

$$
D^{k+2} Q p^{i}=D^{k+2} q^{i}=D^{k+1} q^{i}=D^{k+1} Q p^{i}=p^{i}, \quad 0 \leq i \leq k
$$

Taking into account also Eq. (7.17), this proves Eq. (7.16) for $k+1$. Q.E.D.
It is also interesting to note that the sequence $\left\{x^{k}\right\}$ in Prop. 1.7.2 is identical to the one that would be generated by the preconditioned conjugate gradient method with scaling matrix $H=D^{0}$; i.e., for $k=0,1, \ldots, n-1$, the vector $x^{k+1}$ minimizes $f$ over the linear manifold

$$
M^{k}=\left\{z \mid z=x^{0}+\gamma^{0} D^{0} \nabla f\left(x^{0}\right)+\cdots+\gamma^{k} D^{0} \nabla f\left(x^{k}\right), \gamma^{0}, \ldots, \gamma^{k} \in \Re\right\}
$$

---

This can be proved for the case where $D^{0}=I$ by verifying through induction that for all $k$ there exist scalars $\beta_{i j}^{k}$ such that

$$
D^{k}=I+\sum_{i=0}^{k} \sum_{j=0}^{k} \beta_{i j}^{k} \nabla f\left(x^{i}\right) \nabla f\left(x^{j}\right)^{\prime}
$$

Therefore, for some scalars $b_{i}^{k}$ and all $k$, we have

$$
d^{k}=-D^{k} \nabla f\left(x^{k}\right)=\sum_{i=0}^{k} b_{i}^{k} \nabla f\left(x^{i}\right)
$$

Hence, for all $i, x^{i+1}$ lies on the manifold

$$
M^{i}=\left\{z \mid z=x^{0}+\gamma^{0} \nabla f\left(x^{0}\right)+\cdots+\gamma^{i} \nabla f\left(x^{i}\right), \gamma^{0}, \ldots, \gamma^{i} \in \Re\right\}
$$

and since, by Prop. 1.7.2. the algorithm is a conjugate direction method, $x^{i+1}$ minimizes $f$ over $M^{i}$ based on the results of the preceding section [cf. Eqs. (5.7) and (5.8)]. Thus, when $D^{0}=I$, the algorithm satisfies the defining property of the conjugate gradient method (for all $i, x^{i+1}$ is the unique minimum of $f$ over $M^{i}$ ).

For the case where $D^{0} \neq I$, the proof follows by making a transformation of variables so that in the transformed space the initial matrix is the identity. A consequence of this result is that if line minimization is used and the cost is quadratic, the generated iterates do not depend on the values of the scalar $\xi^{k}$. It turns out that this is also true even when the cost is nonquadratic ([Dix72a], [Dix72b]), which is a rather surprising result. Thus the choice of $\xi^{k}$ makes a difference only if the line minimization is inaccurate.

Finally, we note that multiplying the initial matrix $D^{0}$ by a positive scaling factor can have a significant beneficial effect on the behavior of the algorithm in the initial iterations of an $n$-iteration cycle, and also more generally in the case of a nonquadratic problem. A popular choice is to compute

$$
\tilde{D}^{0}=\frac{p^{0^{\prime}} q^{0}}{q^{0^{\prime}} D^{0} q^{0}} D^{0}
$$

once the vector $x^{1}$ (and hence also $p^{0}$ and $q^{0}$ ) has been obtained, and use $\tilde{D}^{0}$ in place of $D^{0}$ in computing $D^{1}$. Sometimes it is beneficial to scale $D^{k}$ even after the first iteration by multiplication with $p^{k^{\prime}} q^{k} / q^{k^{\prime}} D^{k} q^{k}$; see [OrL74], [Ore73], where it is shown that such scaling can improve the convergence rate.

---

# Comparison of Quasi-Newton Methods with Other Methods 

Let us now consider a nonquadratic problem, and compare the QuasiNewton method of Eqs. (7.1)-(7.2), (7.6)-(7.9) with the conjugate gradient method. One advantage of the quasi-Newton method is that when line search is accurate, the algorithm not only tends to generate conjugate directions but also constructs an approximation to the inverse Hessian matrix. As a result, near convergence to a local minimum with positive definite Hessian, it tends to approximate Newton's method thereby attaining a fast convergence rate. It is significant that this property does not depend on the starting matrix $D^{0}$, and as a result it is not usually necessary to periodically restart the method with a steepest descent-type step, which is something that is essential for the conjugate gradient method.

A second advantage is that the quasi-Newton method is not as sensitive to accuracy in the line search as the conjugate gradient method. This has been verified by extensive computational experience and can be substantiated to some extent by analysis. A partial explanation is that, under essentially no restriction on the line search accuracy, the method generates positive definite matrices $D^{k}$ and hence directions of descent (Prop. 1.7.1).

To compare further the conjugate gradient method and the quasiNewton method, we consider their computational requirements per iteration when $n$ is large. The $k$ th iteration of the conjugate gradient method requires computation of the cost function and its gradient (perhaps several times in the course of the line minimization) together with $O(n)$ operations to compute the conjugate direction $d^{k}$ and the next point $x^{k+1}$. The quasi-Newton method requires roughly the same amount of computation for function and gradient evaluations together with $O\left(n^{2}\right)$ operations to compute the matrix $D^{k}$ and the next point $x^{k+1}$. If the computation needed for a function and gradient evaluation is larger or comparable to $O\left(n^{2}\right)$ operations, the quasi-Newton method requires only slightly more computation per iteration than the conjugate gradient method and holds the edge in view of its other advantages mentioned earlier. In problems where a function and gradient evaluation requires computation time much less than $O\left(n^{2}\right)$ operations, the conjugate gradient method is typically preferable. As an example, we will see in Section 1.9, that in optimal control problems where typically $n$ is very large, a function and a gradient evaluation typically requires $O(n)$ operations. For this reason the conjugate gradient method is typically preferable for these problems.

In general, both the conjugate gradient method and the quasi-Newton algorithm require less computation per iteration than Newton's method, which requires a function, gradient, and Hessian evaluation, as well as $O\left(n^{3}\right)$ operations at each step for computing the Newton direction. This is counterbalanced by the faster speed of convergence of Newton's method. Furthermore, in some cases, special structure can be exploited to compute the Newton direction efficiently. For example in optimal control problems,

---

Newton's method typically requires $O(n)$ operations per iteration versus $O\left(n^{2}\right)$ operations for the quasi-Newton method (see Section 1.9).

# E XERCISES 

### 7.1 (Rank One Quasi-Newton Methods)

Suppose that $D^{k}$ is updated according to the formula

$$
D^{k+1}=D^{k}+\frac{\left(p^{k}-D^{k} q^{k}\right) q^{k^{\prime}}}{q^{k^{\prime}} y^{k}}
$$

where $y^{k}$ is any vector such that $q^{k^{\prime}} y^{k} \neq 0$. Show that we have

$$
D^{k+1} q^{k}=p^{k}
$$

Conclude that for a positive definite quadratic problem, after $n$ steps for which $n$ linearly independent increments $q^{0}, \ldots, q^{n-1}$ are obtained, $D^{n}$ is equal to the inverse Hessian.

### 7.2 (Limited Memory BFGS Method [Noc80])

A major drawback of Quasi-Newton methods for large problems is the large storage requirement. This motivates methods that construct the QuasiNewton direction $d^{k}=-D^{k} \nabla f\left(x^{k}\right)$ using only a limited number of the vectors $p^{i}$ and $q^{i}$ (for example, the last $m$ ). This exercise shows one way to do this.
(a) Show that the BFGS updating formula can be written as

$$
D^{k+1}=V^{k^{\prime}} D^{k} V^{k}+\rho^{k} p^{k} p^{k^{\prime}}
$$

where

$$
\rho^{k}=\frac{1}{q^{k} p^{k}} \quad V^{k}=I-\rho^{k} q^{k} \rho^{k^{\prime}}
$$

(b) Show how to calculate the direction $d^{k}=-D^{k} \nabla f\left(x^{k}\right)$ using $D^{0}$ and the past vectors $p^{i}, q^{i}, i=0,1, \ldots, k-1$.

---

# 1.8 NONDERIVATIVE METHODS 

All the gradient methods examined so far require calculation of at least the gradient $\nabla f\left(x^{k}\right)$ and possibly the Hessian matrix $\nabla^{2} f\left(x^{k}\right)$ at each generated point $x^{k}$. In many problems, either these derivatives are not available in explicit form or they are given by very complicated expressions. In such cases, it may be preferable to use the same algorithms as earlier with all unavailable derivatives approximated by finite differences.

First derivatives may be approximated by the forward difference formula

$$
\frac{\partial f\left(x^{k}\right)}{\partial x^{i}} \approx \frac{1}{h}\left(f\left(x^{k}+h e_{i}\right)-f\left(x^{k}\right)\right)
$$

or by the central difference formula

$$
\frac{\partial f\left(x^{k}\right)}{\partial x^{i}} \approx \frac{1}{2 h}\left(f\left(x^{k}+h e_{i}\right)-f\left(x^{k}-h e_{i}\right)\right)
$$

In these relations, $h$ is a small positive scalar and $e_{i}$ is the $i$ th unit vector (ith column of the identity matrix). In some cases the same value of $h$ can be used for all partial derivatives, but in other cases, particularly when the problem is poorly scaled, it is essential to use a different value of $h$ for each partial derivative. This is a tricky process that often requires trial and error; see the following discussion.

The central difference formula requires twice as much computation as the forward difference formula. However, it is much more accurate. This can be seen by forming the corresponding Taylor series expansions, and by verifying that (in exact arithmetic) the absolute value of the error between the approximation and the actual derivatives is $O(h)$ for the forward difference formula, while it is $O\left(h^{2}\right)$ for the central difference formula. Note that if the central difference formula is used, one obtains at essentially no extra cost an approximation of each diagonal element of the Hessian using the formula

$$
\frac{\partial^{2} f\left(x^{k}\right)}{\left(\partial x_{i}\right)^{2}} \approx \frac{1}{h^{2}}\left(f\left(x^{k}+h e_{i}\right)+f\left(x^{k}-h e_{i}\right)-2 f\left(x^{k}\right)\right)
$$

These approximations can be used in schemes based on diagonal scaling.
To reduce the approximation error, we would like to choose the finite difference interval $h$ as small as possible. Unfortunately, there is a limit on how much $h$ can be reduced due to the roundoff error that occurs when quantities of similar magnitude are subtracted by the computer. In particular, an error $\delta$ due to finite precision arithmetic in evaluating the numerator in Eq. (8.1) [or Eq. (8.2)], results in an error of $\delta / h$ (or $\delta / 2 h$, respectively) in the first derivative evaluation. Roundoff error is particularly evident in the approximate formulas (8.1) and (8.2) near a stationary

---

point where $\nabla f$ is nearly zero, and the relative error size in the gradient approximation becomes very large.

Practical experience suggests that a good policy is to keep the scalar $h$ for each derivative at a fixed value, which roughly balances the approximation error against the roundoff error. Based on the preceding calculations, this leads to the guideline

$$
\begin{aligned}
& \frac{\delta}{h}=O(h) \quad \text { or } \quad h=O\left(\delta^{1 / 2}\right), \quad \text { for the forward difference formula (8.1), } \\
& \frac{\delta}{2 h}=O\left(h^{2}\right) \quad \text { or } \quad h=O\left(\delta^{1 / 3}\right), \quad \text { for the central difference formula (8.2), }
\end{aligned}
$$

where $\delta$ is the error due to finite precision arithmetic in evaluating the numerator in Eq. (8.1) [or Eq. (8.2)]. Thus, a much larger value of $h$ can be used in conjunction with the central difference formula. A good practical rule is to use the forward formula (8.1) until the absolute value of the corresponding approximate derivative becomes less than a certain tolerance; i.e.,

$$
\left|\frac{f\left(x^{k}+h e_{i}\right)-f\left(x^{k}\right)}{h}\right| \leq \epsilon
$$

where $\epsilon>0$ is some prespecified scalar. At that point, a switch to the central difference formula should be made.

Second derivatives may be approximated by the forward difference formula

$$
\frac{\partial^{2} f\left(x^{k}\right)}{\partial x_{i} \partial x_{j}} \approx \frac{1}{h}\left(\frac{\partial f\left(x^{k}+h e_{j}\right)}{\partial x^{i}}-\frac{\partial f\left(x^{k}\right)}{\partial x_{i}}\right)
$$

or the central difference formula

$$
\frac{\partial^{2} f\left(x^{k}\right)}{\partial x_{i} \partial x_{j}} \approx \frac{1}{2 h}\left(\frac{\partial f\left(x^{k}+h e_{j}\right)}{\partial x_{i}}-\frac{\partial f\left(x^{k}-h e_{j}\right)}{\partial x_{i}}\right)
$$

Practical experience suggests that in discretized forms of Newton's method, extreme accuracy in approximating second derivatives is not very important in terms of rate of convergence. For this reason, exclusive use of the forward difference formula (8.3) is adequate in most cases. However, one should certainly check for positive definiteness of the discretized Hessian approximation and introduce modifications if necessary, as discussed in Section 1.4.

# 1.8.1 Coordinate Descent 

There are several other nonderivative methods for minimizing differentiable functions. A particularly important algorithm is the coordinate descent method. Here the cost is minimized along one coordinate direction at each iteration. The order in which coordinates are chosen may vary in

---

the course of the algorithm. In the case where this order is cyclical, given $x^{k}$, the $i$ th coordinate of $x^{k+1}$ is determined by

$$
x_{i}^{k+1}=\arg \min _{\xi \in \Re} f\left(x_{1}^{k+1}, \ldots, x_{i-1}^{k+1}, \xi, x_{i+1}^{k}, \ldots, x_{n}^{k}\right)
$$

see Fig. 1.8.1. The method can also be used for minimization of $f$ subject to upper and lower bounds on the variables $x^{i}$ (the minimization over $\xi \in \Re$ in the preceding equation is replaced by minimization over the appropriate interval). We will analyze the method within this more general context in the next chapter.

An important advantage of the coordinate descent method is that it is well suited for parallel computation. In particular, suppose that there is a subset of coordinates $x_{i_{1}}, x_{i_{2}}, \ldots, x_{i_{m}}$, which are not coupled through the cost function, that is, $f(x)$ can be written as $\sum_{r=1}^{m} f_{i_{r}}(x)$, where for each $r, f_{i_{r}}(x)$ does not depend on the coordinates $x_{i_{s}}$ for all $s \neq r$. Then one can perform the $m$ coordinate descent iterations

$$
x_{i_{r}}^{k+1}=\arg \min _{\xi} f_{i_{r}}\left(x^{k}+\xi e_{i_{r}}\right), \quad r=1, \ldots, m
$$

independently and in parallel. Thus, in problems with special structure where the set of coordinates can be partitioned into $p$ subsets with the independence property just described, one can perform a full cycle of coordinate descent iterations in $p$ (as opposed to $n$ ) parallel steps (assuming of course that a sufficient number of parallel processors is available).

The coordinate descent method generally has similar convergence properties to steepest descent. For continuously differentiable cost functions, it can be shown to generate sequences whose limit points are stationary, although the proof of this is sometimes complicated and requires some additional assumptions (see Prop. 2.7.1 in Section 2.7, which deals with a constrained version of coordinate descent and requires strict convexity of the cost function along each coordinate). There is also a great deal of analysis of coordinate descent in a context where its use is particularly favorable, namely in solving dual problems (see Section 6.2). Within this context, the strict convexity assumption is neither satisfied nor is it essential (see the references given in Chapter 6). The convergence rate of coordinate descent to nonsingular and singular local minima can be shown to be linear and sublinear, respectively, similar to steepest descent. Often, the choice between coordinate descent and steepest descent is dictated by the structure of the cost function. Both methods can be very slow, but for many practical contexts, they can be quite effective.

# 1.8.2 Direct Search Methods 

In the coordinate descent method we search along the fixed set of coordinate directions and we are guaranteed a cost improvement at a nonstationary point because these directions are linearly independent. This

---

![[chapter_21_p145_img28.jpeg]]

Figure 1.8.1. Illustration of the coordinate descent method.
idea can be generalized by using a different set of directions and by occasionally changing this set of directions with the aim of accelerating convergence. There are a number of methods of this type: Rosenbrock's method [Ros60a], the pattern search algorithm of Hooke and Jeeves [HoJ61], and the simplex algorithms of Spendley, Hext, and Himsworth [SHH62], and Nelder and Mead [NeM65]. Unfortunately, the rationale of these methods often borders on the heuristic, and their theoretical convergence properties are often unsatisfactory. However, these methods are often fairly simple to implement and like the coordinate descent method, they do not require gradient calculations. We describe the Nelder and Mead simplex method (not to be confused with the simplex method of linear programming), which has enjoyed considerable popularity.

At the typical iteration of this method, we start with a simplex, that is, the convex hull of $n+1$ points, $x^{0}, x^{1}, \ldots, x^{n}$, and we end up with another simplex. Let $x_{\text {min }}$ and $x_{\text {max }}$ denote the "best" and "worst" vertices of the simplex, that is the vertices satisfying

$$
\begin{aligned}
& f\left(x_{\min }\right)=\min _{i=0,1, \ldots, n} f\left(x^{i}\right) \\
& f\left(x_{\max }\right)=\max _{i=0,1, \ldots, n} f\left(x^{i}\right)
\end{aligned}
$$

Let also $\hat{x}$ denote the centroid of the face of the simplex formed by the vertices other than $x_{\max }$

$$
\hat{x}=\frac{1}{n}\left(-x_{\max }+\sum_{i=0}^{n} x^{i}\right)
$$

The iteration replaces the worst vertex $x_{\max }$ by a "better" one. In particular, the reflection point $x_{\text {ref }}=2 x-x_{\text {max }}$ is computed, which lies on the line passing through $x_{\text {max }}$ and $x$, and is symmetric to $x_{\text {max }}$ with respect to $\hat{x}$. Depending on the cost value of $x_{\text {ref }}$ relative to the points of the simplex other than $x_{\text {max }}$, a new vertex $x_{\text {new }}$ is computed, and a new simplex is formed from the old by replacing the vertex $x_{\text {max }}$ by $x_{\text {new }}$, while keeping the other $n$ vertices.

---

# Typical Iteration of the Simplex Method 

Step 1: (Reflection Step) Compute

$$
x_{r e f}=2 \hat{x}-x_{\max }
$$

Then compute $x_{\text {new }}$ according to the following three cases:
(1) ( $x_{\text {ref }}$ has min cost) If $f\left(x_{\text {min }}\right)>f\left(x_{\text {ref }}\right)$, go to Step 2.
(2) ( $x_{\text {ref }}$ has intermediate cost) If $\max \left\{f\left(x^{i}\right) \mid x^{i} \neq x_{\max }\right\}>f\left(x_{\text {ref }}\right) \geq$ $f\left(x_{\text {min }}\right)$, go to Step 3.
(3) ( $x_{\text {ref }}$ has max cost) If $f\left(x_{\text {ref }}\right) \geq \max \left\{f\left(x^{i}\right) \mid x^{i} \neq x_{\max }\right\}$, go to Step 4.

Step 2: (Attempt Expansion) Compute

$$
x_{\exp }=2 x_{r e f}-\hat{x}
$$

Define

$$
x_{\text {new }}=\left\{\begin{array}{ll}
x_{\exp } & \text { if } f\left(x_{\exp }\right)<f\left(x_{\text {ref }}\right) \\
x_{\text {ref }} & \text { otherwise }
\end{array}\right.
$$

and form the new simplex by replacing the vertex $x_{\text {max }}$ with $x_{\text {new }}$.
Step 3: (Use Reflection) Define $x_{\text {new }}=x_{\text {ref }}$, and form the new simplex by replacing the vertex $x_{\text {max }}$ with $x_{\text {new }}$.

Step 4: (Perform Contraction) Define

$$
x_{\text {new }}=\left\{\begin{array}{ll}
\frac{1}{2}\left(x_{\max }+\hat{x}\right) & \text { if } f\left(x_{\max }\right) \leq f\left(x_{\text {ref }}\right) \\
\frac{1}{2}\left(x_{\text {ref }}+\hat{x}\right) & \text { otherwise }
\end{array}\right.
$$

and form the new simplex by replacing the vertex $x_{\text {max }}$ with $x_{\text {new }}$.
The reflection step and the subsequent possible steps of the iteration are illustrated in Fig. 1.8.2 (a)-(d). The entire method is illustrated in Fig. 1.8.3. Exercise 8.3 shows a cost improvement property of the method in the case where $f$ is strictly convex. However, there are no known convergence results for the method. Furthermore, when the cost function is not convex, it is possible that the new simplex vertex $x_{\text {new }}$ has larger cost value than the old vertex $x_{\text {max }}$. In this case a modification that has been suggested is to "shrink" the old simplex towards the best vertex $x_{\text {min }}$, that is, form a new simplex by replacing all the vertices $x^{i}, i=0,1, \ldots, n$, by

$$
\bar{x}^{i}=\frac{1}{2}\left(x^{i}+x_{\min }\right), \quad i=0,1, \ldots, n
$$

The method as given above seems to work reasonably well in practice, particularly for problems of relatively small dimension (say up to 10). However, it is not guaranteed to have desirable convergence properties, and in fact a convergence counterexample is given in [McK94]. Reference

---

[Tse95a] provides a relatively simple modification with satisfactory convergence properties. There are also a number of related methods, some of which have demonstrable convergence properties; see [DeT91] and [Tor91]. Note that the constants used in Eqs. (8.9), (8.10), and (8.11) are somewhat arbitrary, as suggested by the interpretation of the method given in Figs. 1.8.2 and 1.8.3. More general forms of these equations are

$$
\begin{gathered}
x_{r e f}=\hat{x}+\beta\left(\hat{x}-x_{\max }\right) \\
x_{\text {exp }}=x_{r e f}+\gamma\left(x_{r e f}-\hat{x}\right) \\
x_{c o n}= \begin{cases}\theta x_{\max }+(1-\theta) \hat{x} & \text { if } f\left(x_{\max }\right) \leq f\left(x_{r e f}\right) \\
\theta x_{r e f}+(1-\theta) \hat{x} & \text { otherwise }\end{cases}
\end{gathered}
$$

where $\beta>0, \gamma>0$, and $\theta \in(0,1)$ are scalars known as the reflection coefficient, the expansion coefficient, and the contraction coefficient, respectively. The formulas of Eqs. (8.9)-(8.11) correspond to $\beta=1, \gamma=1$, and $\theta=1 / 2$, respectively.
![[chapter_21_p147_img29.jpeg]]

Figure 1.8.2. Illustration of the reflection step and the possible subsequent steps of an iteration of the simplex method. In (a), the new vertex $x_{\text {new }}$ is determined via the expansion Step 2. In (b), $x_{\text {new }}$ is determined via Step 3, and the reflection step is accepted. In (c) and (d), $x_{\text {new }}$ is determined via the contraction Step 4.

---

![[chapter_21_p148_img30.jpeg]]

Figure 1.8.3. Illustration of three iterations of the simplex method, which generate the points $x^{3}, x^{4}$, and $x^{5}$, starting from the simplex $x^{0}, x^{1}, x^{2}$. The simplex obtained after these three iterations consists of $x^{2}, x^{4}$, and $x^{5}$.

# E XERCISES 

## 8.1

Let $f: \Re^{n} \mapsto \Re$ be continuously differentiable, let $p_{1}, \ldots, p_{n}$ be linearly independent vectors, and suppose that for some $x^{*}, \alpha=0$ is a stationary point of each of the one-dimensional functions $g_{i}(\alpha)=f\left(x^{*}+\alpha p_{i}\right), i=1, \ldots, n$. Show that $x^{*}$ is a stationary point of $f$.

### 8.2 (Stepsize Selection in Jacobi Methods)

Let $f: \Re^{n} \rightarrow \Re$ be a continuously differentiable convex function. For a given $x \in \Re^{n}$ and all $i=1, \ldots, n$, define the vector $\bar{x}$ by

$$
\bar{x}_{i}=\arg \min _{\xi \in \Re} f\left(x_{1}, \ldots, x_{i-1}, \xi, x_{i+1}, \ldots, x_{n}\right)
$$

The Jacobi method is defined by the iteration

$$
x:=x+\alpha(\bar{x}-x)
$$

where $\alpha$ is a positive stepsize parameter.
(a) Show that if $x$ does not minimize $f$, then the Jacobi iteration reduces the value of $f$ when $\alpha=1 / n$.

---

(b) [Rus95] Consider the case where $f$ has the form

$$
f\left(x_{i}\right)=\sum_{j=1}^{J} f_{j}\left(\sum_{i=1}^{n} a_{i}, x_{i}\right)
$$

where $a_{j i}$ are given scalars. Let $m=\max _{j=1, \ldots, J} m_{j}$, where $m_{j}$ is the number of indices $i$ for which $a_{j i}$ is nonzero. Show that if $x$ does not minimize $f$, the Jacobi iteration reduces the value of $f$ when $\alpha=1 / m$.

# 8.3 

Consider the simplex method applied to a strictly convex function $f$. Show that at each iteration, either $f\left(x_{\text {max }}\right)$ decreases strictly, or else the number of vertices $x^{\prime}$ of the simplex such that $f\left(x^{\prime \prime}\right)=f\left(x_{\text {max }}\right)$ decreases by at least one.

### 1.9 DISCRETE-TIME OPTIMAL CONTROL PROBLEMS*

In this section, we consider a class of optimization problems involving a discrete-time dynamic system, that is, a vector difference equation. Such problems arise often in applications and are frequently challenging because of their large dimension. Continuous-time versions of these problems are the focus of the modern theory of the calculus of variations and the Pontryagin maximum principle, and give rise to some of the most fascinating mathematical problems of optimization. We will focus on some structural aspects of optimal control problems, which can be effectively exploited both in analysis and computation.

Let us consider the problem of finding sequences $u=\left(u_{0}, u_{1}, \ldots, u_{N-1}\right)$ and $x=\left(x_{1}, x_{2}, \ldots, x_{N}\right)$, which minimize

$$
g_{N}\left(x_{N}\right)+\sum_{i=0}^{N-1} g_{i}\left(x_{i}, u_{i}\right)
$$

subject to the constraints

$$
\begin{gathered}
x_{i+1}=f_{i}\left(x_{i}, u_{i}\right), \quad i=0, \ldots, N-1 \\
x_{0}: \text { given, } \\
x_{i} \in X_{i} \subset \Re^{n}, \quad i=1, \ldots, N \\
u_{i} \in U_{i} \subset \Re^{m}, \quad i=0, \ldots, N-1
\end{gathered}
$$

---

We refer to $u_{i}$ as the control vectors, and to $u=\left(u_{0}, u_{1}, \ldots, u_{N-1}\right)$ as a control trajectory. We refer to $x_{i}$ as the state vectors, and to the sequence $x=\left(x_{0}, x_{1}, \ldots, x_{N}\right)$ as a state trajectory. The equation $x_{i+1}=f_{i}\left(x_{i}, u_{i}\right)$ is called the system equation and specifies uniquely the state trajectory, which corresponds to a given control trajectory. The functions $f_{i}: \Re^{n} \times \Re^{m} \mapsto \Re^{n}$ are given and they will be assumed once or twice differentiable for the most part of the following. The sets $X_{i}$ and $U_{i}$ specify constraints on the state and control vectors, respectively. They are usually represented by equations or inequalities. However, for the time being we will not specify them further. A more general constraint, that we will discuss on occasion has the form

$$
\left(x_{i}, u_{i}\right) \in \Omega_{i} \subset \Re^{n} \times \Re^{m}, \quad i=0, \ldots, N-1
$$

The functions $g_{N}: \Re^{n} \mapsto \Re$ and $g_{i}: \Re^{n} \times \Re^{m} \mapsto \Re, i=0, \ldots, N-1$, are given and for the most part, they will be assumed to be once or twice differentiable.

# Example 9.1 (Reservoir Regulation) 

Let $x_{i}$ denote the volume of water held in a reservoir at the $i$ th of $N$ time periods. The volume $x_{i}$ evolves according to

$$
x_{i+1}=x_{i}-u_{i}, \quad \forall i=0, \ldots, N-1
$$

where $u_{i}$ is water used for some productive purpose in period $i$. This is the system equation, with the volume $x_{i}$ viewed as the state and the outflow $u_{i}$ viewed as the control. There is a cost $G\left(x_{N}\right)$ for the terminal state being $x_{N}$ and there is a cost $g_{i}\left(u_{i}\right)$ for outflow $u_{i}$ at period $i$. For example, when $u_{i}$ is used for electric power generation, $g_{i}\left(u_{i}\right)$ may be equal to minus the value of power produced from $u_{i}$. We want to choose the outflows $u_{0}, u_{1}, \ldots, u_{N-1}$ so as to minimize

$$
G\left(x_{N}\right)+\sum_{i=0}^{N-1} g_{i}\left(u_{i}\right)
$$

while observing some constraints on the volume $x_{i}$ (e.g. $x_{i}$ should lie between some upper and lower bounds) and some constraints on the outflow $u_{i}$ (e.g. $u_{i} \geq 0$ ).

There are also multidimensional versions of the problem, involving several reservoirs that are interconnected in the sense that the outflow $u_{i}$ from one becomes inflow to another (in addition to serving some other productive purpose). This leads to an optimal control formulation where the state and control vectors are multidimensional with dimension equal to the number of reservoirs.

Let us consider the case where there are no state constraints, i.e.,

$$
X_{i}=\Re^{n}, \quad i=1, \ldots, N
$$

---

Then, one may reduce problem (9.1)-(9.4) (which is a constrained problem in $x$ and $u$ ) to a problem which involves only the control variables $u_{0}, \ldots, u_{N-1}$. To see this, note that to any given control trajectory

$$
u=\left(u_{0}, u_{1}, \ldots, u_{N-1}\right)
$$

there corresponds a unique state trajectory via the system equation $x_{i+1}=$ $f_{i}\left(x_{i}, u_{i}\right)$. We may write this correspondence abstractly as

$$
x_{i}=\phi_{i}(u)=\phi_{i}\left(u_{0}, \ldots, u_{N-1}\right), \quad i=1, \ldots, N
$$

where $\phi_{i}$ are appropriate functions, determined by the functions $f_{i}$. (Actually in each of the functions $\phi_{i}$, only the variables $u_{0}, \ldots, u_{i-1}$ enter explicitly. However, our notation is technically correct and will prove convenient.)

We may substitute $x_{i}$, using the functions $\phi_{i}$ of Eq. (9.5), into the cost function (9.1) and write the problem as

$$
\begin{aligned}
& \text { minimize } J(u)=g_{N}\left(\phi_{N}(u)\right)+\sum_{i=0}^{N-1} g_{i}\left(\phi_{i}(u), u_{i}\right) \\
& \text { subject to } u_{i} \in U_{i}, \quad i=0, \ldots, N-1
\end{aligned}
$$

If the controls are also unconstrained, that is,

$$
U_{i}=\Re^{m}, \quad i=0, \ldots, N-1
$$

problem (9.6) is an unconstrained minimization problem of the type that we have examined in this chapter. Thus, once we calculate the gradient and Hessian matrix of $J$, we can write explicitly necessary conditions and sufficient conditions for optimality.

We will calculate the derivatives of $J$ in two different ways, both of which are valuable, since they provide complementary insights.

# Calculation of $\nabla J(u)$ (1st Method) 

Let us first assume that the cost function has the form

$$
J(u)=g_{N}\left(x_{N}\right)=g_{N}\left(\phi_{N}(u)\right)
$$

i.e., there is only a terminal state cost. As we shall see shortly, the general case can be reduced to this case.

We have for every $i=0, \ldots, N-1$, using the chain rule,

$$
\nabla_{u_{i}} J(u)=\nabla_{u_{i}} \phi_{N} \cdot \nabla g_{N}
$$

---

where the derivatives are calculated along the current control trajectory $u$ and the corresponding state trajectory $x$, and where $x$ and $u$ are related by

$$
x_{i}=\phi_{i}(u), \quad i=1, \ldots, N
$$

[cf. Eq. (9.5)]. From Eq. (9.7), we obtain, using again the chain rule,

$$
\begin{aligned}
\nabla_{u_{i}} J(u) & =\nabla_{u_{i}} x_{i+1} \cdot \nabla_{x_{i+1}} x_{i+1} \cdots \nabla_{x_{N-1}} x_{N} \cdot \nabla g_{N} \\
& =\nabla_{u_{i}} f_{i} \cdot \nabla_{x_{i+1}} f_{i+1} \cdots \nabla_{x_{N-1}} f_{N-1} \cdot \nabla g_{N}
\end{aligned}
$$

By defining the, so called, costate vectors $p_{i} \in \Re^{n}, i=1, \ldots, N$, via the equations

$$
\begin{gathered}
p_{i}=\nabla_{x_{i}} f_{i} \cdots \nabla_{x_{N-1}} f_{N-1} \cdot \nabla g_{N}, \quad i=1, \ldots, N-1 \\
p_{N}=\nabla g_{N}
\end{gathered}
$$

we obtain

$$
\nabla_{u_{i}} J(u)=\nabla_{u_{i}} f_{i} \cdot p_{i+1}, \quad i=0, \ldots, N-1
$$

Note also that $p_{1}, \ldots, p_{N}$ are generated (backwards in time) by means of the following equation, known as the adjoint equation,

$$
p_{i}=\nabla_{x_{i}} f_{i} \cdot p_{i+1}, \quad i=1, \ldots, N-1
$$

starting from

$$
p_{N}=\nabla g_{N}
$$

Consider now the general case of the cost function

$$
J(u)=g_{N}\left(x_{N}\right)+\sum_{i=0}^{N-1} g_{i}\left(x_{i}, u_{i}\right)
$$

which involves the intermediate cost terms $g_{i}\left(x_{i}, u_{i}\right)$. We may reduce the corresponding problem to one which involves a terminal state cost only, by introducing an additional state variable $y_{i}$ and a corresponding equation by

$$
\begin{gathered}
y_{i+1}=y_{i}+g_{i}\left(x_{i}, u_{i}\right) \\
y_{0}=0
\end{gathered}
$$

Defining a new state vector

$$
\tilde{x}_{i}=\binom{x_{i}}{y_{i}}
$$

and a new system equation

$$
\tilde{x}_{i+1}=\binom{f_{i}\left(x_{i}, u_{i}\right)}{y_{i}+g_{i}\left(x_{i}, u_{i}\right)} \equiv \tilde{f}_{i}\left(\tilde{x}_{i}, u_{i}\right)
$$

---

the cost function (9.13) takes the form

$$
J(u)=g_{N}\left(x_{N}\right)+y_{N} \equiv \tilde{g}_{N}\left(\tilde{x}_{N}\right)
$$

This cost function involves only a terminal state cost and is of the type examined previously. We have from Eqs. (9.10)-(9.12),

$$
\nabla_{u_{i}} J(u)=\nabla_{u_{i}} \tilde{f}_{i} \cdot \tilde{p}_{i+1}, \quad i=0, \ldots, N-1
$$

where

$$
\begin{gathered}
\tilde{p}_{i}=\nabla_{\tilde{x}_{i}} \tilde{f}_{i} \cdot \tilde{p}_{i+1}, \quad i=1, \ldots, N-1 \\
\tilde{p}_{N}=\nabla \tilde{g}_{N}
\end{gathered}
$$

From Eqs. (9.14)-(9.16), we have by writing

$$
\tilde{p}_{i}=\binom{p_{i}}{z_{i}}, \quad p_{i} \in \Re^{n}, z_{i} \in \Re
$$

the following form of the adjoint equation

$$
\begin{aligned}
\tilde{p}_{N} & =\binom{p_{N}}{z_{N}}=\binom{\nabla g_{N}}{1} \\
\tilde{p}_{i}=\binom{p_{i}}{z_{i}} & =\left(\begin{array}{cc}
\nabla_{x_{i}} f_{i} & \nabla_{x_{i}} g_{i} \\
0 & 1
\end{array}\right)\binom{p_{i+1}}{z_{i+1}}
\end{aligned}
$$

As a result, we have

$$
z_{i}=1, \quad i=1, \ldots, N
$$

and we obtain the final form of the adjoint equation

$$
\begin{gathered}
p_{i}=\nabla_{x_{i}} f_{i} \cdot p_{i+1}+\nabla_{x_{i}} g_{i} \\
p_{N}=\nabla g_{N}
\end{gathered}
$$

Furthermore, from Eq. (9.17),

$$
\nabla_{u_{i}} J(u)=\left(\nabla_{u_{i}} f_{i} \quad \nabla_{u_{i}} g_{i}\right)\binom{p_{i+1}}{z_{i+1}}
$$

and finally

$$
\nabla_{u_{i}} J(u)=\nabla_{u_{i}} f_{i} \cdot p_{i+1}+\nabla_{u_{i}} g_{i}
$$

Equation (9.20) yields the cost gradient in terms of the costate vectors obtained from the adjoint equation (9.18) and (9.19), with all derivatives evaluated along the control trajectory $u$ and the corresponding state trajectory under consideration. Thus to compute the gradient $\nabla J(u)$ :
(1) We calculate the state trajectory corresponding to $u$ by forwards propagation of the system equation $x_{i+1}=f_{i}\left(x_{i}, u_{i}\right)$.
(2) We generate the costate vectors by backwards propagation of the adjoint equation Eqs. (9.18), starting from the terminal condition $p_{N}=\nabla g_{N}$. (All partial derivatives are evaluated at the current state and control trajectories.)
(3) We use Eq. (9.20) to compute the components $\nabla_{u_{i}} J(u)$ of the gradient $\nabla J(u)$.

---

# Calculation of $\nabla J(u)$ and $\nabla^{2} J(u)$ (2nd Method) 

The alternative method for calculating the derivatives of $J$ is based on writing all the system equations $x_{i+1}=f_{i}\left(x_{i}, u_{i}\right)$ compactly as

$$
h(\phi(u), u)=0
$$

Here, $u=\left(u_{0}, u_{1}, \ldots, u_{N-1}\right)$ is the control trajectory and $\phi(u)$ is the corresponding state trajectory

$$
x=\phi(u)
$$

where the function $\phi$ maps $\Re^{N m}$ into $\Re^{N n}$ and is defined in terms of the functions $\phi_{i}$ of Eq. (9.5) as

$$
\phi(u)=\left(\begin{array}{c}
\phi_{1}(u) \\
\vdots \\
\phi_{N}(u)
\end{array}\right)=\left(\begin{array}{c}
x_{1} \\
\vdots \\
x_{N}
\end{array}\right)
$$

The equation $h(\phi(u), u)=0$ is a compact representation of the equations

$$
\begin{gathered}
f_{0}\left(x_{0}, u_{0}\right)-x_{1}=0 \\
f_{1}\left(x_{1}, u_{1}\right)-x_{2}=0 \\
\vdots \\
f_{N-1}\left(x_{N-1}, u_{N-1}\right)-x_{N}=0
\end{gathered}
$$

and is satisfied by every control trajectory $u$.
Similarly the cost function may be written abstractly as

$$
F(\phi(u), u)=g_{N}\left(\phi_{N}(u)\right)+\sum_{i=0}^{N-1} g_{i}\left(\phi_{i}(u), u_{i}\right)
$$

and the optimal control problem becomes

$$
\begin{aligned}
& \operatorname{minimize} \quad J(u)=F(\phi(u), u) \\
& \text { subject to } u \in \Re^{N m}
\end{aligned}
$$

To obtain the derivatives of $J(u)$, we use the trick of writing for any vector $p \in \Re^{N n}$,

$$
J(u)=F(\phi(u), u)+h(\phi(u), u)^{\prime} p
$$

[cf. Eq. (9.21)] and we select cleverly $p$ so that the derivatives of $J(u)$ are easily calculated. We have

$$
\begin{aligned}
\nabla J(u)= & \nabla \phi(u)\left(\nabla_{x} F(\phi(u), u)+\nabla_{x} h(\phi(u), u) \cdot p\right) \\
& +\nabla_{u} F(\phi(u), u)+\nabla_{u} h(\phi(u), u) \cdot p
\end{aligned}
$$

---

where in the preceding equation, $\nabla_{x}$ denotes gradient with respect to the first argument $[x=\phi(u)], \nabla_{u}$ denotes gradient with respect to the second argument $u$, and $\nabla \phi(u)$ is the $N m \times N n$ gradient matrix of $\phi$ evaluated at $u$, i.e.,

$$
\nabla \phi(u)=\left[\nabla_{u} \phi_{1}(u) \vdots \nabla_{u} \phi_{2}(u) \vdots \cdots \vdots \nabla_{u} \phi_{N}(u)\right]
$$

Equation (9.25) holds for every $p \in \Re^{N n}$. Suppose that $p$ is selected so that it satisfies the equation

$$
\nabla_{x} F(\phi(u), u)+\nabla_{x} h(\phi(u), u) \cdot p(u)=0
$$

where we denote by $p(u)$ the particular value of $p$ satisfying the equation above. We will see shortly that Eq. (9.26) has a unique solution with respect to $p$, so our notation is justified. We have then from Eq. (9.25),

$$
\nabla J(u)=\nabla_{u} F(\phi(u), u)+\nabla_{u} h(\phi(u), u) \cdot p(u)
$$

For notational convenience, we introduce the function $L: \Re^{N n} \times \Re^{N m} \times$ $\Re^{N n} \mapsto \Re$ given by

$$
L(x, u, p)=F(x, u)+h(x, u)^{\prime} p
$$

We then have

$$
\nabla J(u)=\nabla_{u} L(\phi(u), u, p(u))
$$

To compute the Hessian matrix of the cost, we write Eq. (9.25) as

$$
\nabla J(u)=\nabla \phi(u) \nabla_{x} L(\phi(u), u, p)+\nabla_{u} L(\phi(u), u, p)
$$

and we differentiate this equation with respect to $u$. By using also the equation

$$
\nabla_{x} L(\phi(u), u, p(u))=0
$$

[cf. Eq. (9.26)] we obtain for $p=p(u)$

$$
\begin{aligned}
\nabla^{2} J(u)= & \nabla \phi(u) \nabla_{x x}^{2} L(\phi(u), u, p(u)) \nabla \phi(u)^{\prime}+\nabla \phi(u) \nabla_{x u}^{2} L(\phi(u), u, p(u)) \\
& +\nabla_{u x}^{2} L(\phi(u), u, p(u)) \nabla \phi(u)^{\prime}
\end{aligned}
$$

In order to calculate the derivatives of $J$, we must calculate the solution $p(u)$ of Eq. (9.26) and then evaluate the derivatives of $L$. We have:

$$
L(x, u, p)=g_{N}\left(x_{N}\right)+\sum_{i=0}^{N-1} g_{i}\left(x_{i}, u_{i}\right)+\sum_{i=1}^{N}\left(f_{i-1}\left(x_{i-1}, u_{i-1}\right)-x_{i}\right)^{\prime} p_{i}
$$

---

Equation (9.26) can be written as $\nabla_{x} L(\phi(u), u, p(u))=0$, or by using Eq. $(9.31)$,

$$
\begin{aligned}
& \left(\begin{array}{c}
\nabla g_{N} \\
\nabla_{x_{N-1}} g_{N-1} \\
\nabla_{x_{N-2}} g_{N-2} \\
\vdots \\
\nabla_{x_{1}} g_{1}
\end{array}\right) \\
& +\left(\begin{array}{ccccc}
-I & 0 & 0 & \cdots & 0 \\
\nabla_{x_{N-1}} f_{N-1} & -I & 0 & \cdots & 0 \\
0 & \nabla_{x_{N-2}} f_{N-2} & -I & \cdots & 0 \\
\vdots & \vdots & \vdots & \cdots & \vdots \\
0 & 0 & 0 & \cdots & -I
\end{array}\right)\left(\begin{array}{c}
p_{N} \\
p_{N-1} \\
p_{N-2} \\
\vdots \\
p_{1}
\end{array}\right)=0 .
\end{aligned}
$$

By rewriting this equation, we have

$$
\begin{gathered}
p_{i}=\nabla_{x_{i}} f_{i} \cdot p_{i+1}+\nabla_{x_{i}} g_{i}, \quad i=1, \ldots, N-1 \\
p_{N}=\nabla g_{N}
\end{gathered}
$$

which is the same as the adjoint equation derived earlier [Eqs. (9.18) and (9.19)]. Thus, the adjoint equation specifies uniquely $p(u)$.

Also, from Eqs. (9.28), (9.29), and (9.31),

$$
\nabla_{u_{i}} J(u)=\nabla_{u_{i}} f_{i} \cdot p_{i+1}+\nabla_{u_{i}} g_{i}, \quad i=0, \ldots, N-1
$$

which is the same equation for the gradient $\nabla J$ as the one obtained earlier.
Concerning the Hessian matrix $\nabla^{2} J$, it is given by Eqs. (9.30) and (9.31). We only give here the form of the gradient matrix $\nabla \phi(u)$

$$
\nabla \phi(u)=\left(\begin{array}{cccc}
\nabla_{u_{0}} f_{0} & \nabla_{x_{1}} f_{1} \cdot \nabla_{u_{0}} f_{0} & \cdots & \nabla_{x_{N-1}} f_{N-1} \cdots \nabla_{x_{1}} f_{1} \cdot \nabla_{u_{0}} f_{0} \\
0 & \nabla_{u_{1}} f_{1} & \cdots & \nabla_{x_{N-1}} f_{N-1} \cdots \nabla_{x_{2}} f_{2} \cdot \nabla_{u_{1}} f_{1} \\
\vdots & \vdots & \cdots & \vdots \\
0 & 0 & \cdots & \nabla_{u_{N-1}} f_{N-1}
\end{array}\right)
$$

The expression (9.30) for $\nabla^{2} J$ is quite complex, but will be useful when we will discuss Newton's method later in this section.

# Example 9.2 (Neural Networks) 

Recall Example 5.3 in Section 5.1, which describes how the neural network training problem can be posed as a least squares problem. It can be seen that this problem can also be viewed as an optimal control problem, where the

---

weight vector at each stage is the control vector for that stage. Mathematically, the problem is to find the control trajectory $u$ that minimizes

$$
J(u)=\sum_{j=1}^{m} J_{j}(u)
$$

where

$$
J_{j}(u)=\frac{1}{2}\left\|z_{j}-h\left(u, y_{j}\right)\right\|^{2}
$$

$\left(z_{j}, y_{j}\right)$ is the $j$ th input-output pair in the training set, and $h$ is an appropriate function defined by the functional relationship between input and output of the neural network. Each gradient $\nabla J_{j}(u)$ can be calculated as described earlier via the corresponding adjoint equation, and the cost gradient is obtained as

$$
\nabla J(u)=\sum_{j=1}^{m} \nabla J_{j}(u)
$$

The process of calculating the gradient by means of the adjoint equation is sometimes called back-propagation and corresponding gradient-like methods are sometimes called back-propagation methods.

# First Order Necessary Condition 

The preceding expressions for the cost gradient can be used to write the first order necessary optimality condition

$$
\nabla J\left(u^{*}\right)=0
$$

for

$$
u^{*}=\left(u_{0}^{*}, u_{1}^{*}, \ldots, u_{N-1}^{*}\right)
$$

to be a local minimum of $J$. It is customary to write this condition in terms of the Hamiltonian function, defined for each $i$ by

$$
H_{i}\left(x_{i}, u_{i}, p_{i+1}\right)=g_{i}\left(x_{i}, u_{i}\right)+p_{i+1}^{\prime} f_{i}\left(x_{i}, u_{i}\right)
$$

By using the expression of $\nabla J\left(u^{*}\right)$ [cf. Eq. (9.32)-(9.34)] we obtain:

---

Proposition 1.9.1: Let $u^{*}=\left(u_{0}^{*}, u_{1}^{*}, \ldots, u_{N-1}^{*}\right)$ be a local minimum control trajectory and let $x^{*}=\left(x_{1}^{*}, x_{2}^{*}, \ldots, x_{N}^{*}\right)$ be the corresponding state trajectory. Then we have

$$
\nabla_{u_{i}} H_{i}\left(x_{i}^{*}, u_{i}^{*}, p_{i+1}^{*}\right)=0, \quad i=0, \ldots, N-1
$$

where the costate vectors $p_{1}^{*}, \ldots, p_{N}^{*}$ are obtained from the adjoint equation

$$
p_{i}^{*}=\nabla_{x_{i}} H_{i}\left(x_{i}^{*}, u_{i}^{*}, p_{i+1}^{*}\right), \quad i=1, \ldots, N-1
$$

with the terminal condition

$$
p_{N}^{*}=\nabla g_{N}\left(x_{N}^{*}\right)
$$

# Example 9.3 (Linear System and Quadratic Cost) 

Consider the case where the system is linear

$$
x_{i+1}=A_{i} x_{i}+B_{i} u_{i}, \quad i=0, \ldots, N-1
$$

and the cost function is quadratic of the form

$$
J(u)=\frac{1}{2}\left\{x_{N}^{\prime} Q_{N} x_{N}+\sum_{i=0}^{N-1}\left(x_{i}^{\prime} Q_{i} x_{i}+u_{i}^{\prime} R_{i} u_{i}\right)\right\}
$$

The $n \times n$ matrices $A_{i}$ and the $n \times m$ matrices $B_{i}$ are given. The matrices $Q_{i}$ are assumed symmetric and positive semidefinite, and the matrices $R_{i}$ are assumed symmetric and positive definite. There are no constraints on the state or control vectors.

Let $u^{*}=\left(u_{0}^{*}, u_{1}^{*}, \ldots, u_{N-1}^{*}\right)$ be an optimal control trajectory and let $x^{*}=\left(x_{1}^{*}, x_{2}^{*}, \ldots, x_{N}^{*}\right)$ be the corresponding state trajectory. The adjoint equation is given by

$$
\begin{gathered}
p_{i}^{*}=A_{i}^{\prime} p_{i+1}^{*}+Q_{i} x_{i}^{*}, \quad i=1, \ldots, N-1 \\
p_{N}^{*}=Q_{N} x_{N}^{*}
\end{gathered}
$$

It can be seen from Eq. (9.39), that the state trajectory $x$ is linearly related to the control trajectory $u$, so the cost function $J(u)$ of Eq. (9.40) is a quadratic function of $u$. Because the matrices $Q_{i}$ and $R_{i}$ are assumed positive semidefinte and positive definite, respectively, $J(u)$ is a convex, positive

---

definite, quadratic function. Therefore, the necessary optimality condition $\nabla J\left(u^{*}\right)=0$ of Prop. 1.8.1 is also sufficient. It can be written as

$$
\nabla_{u_{i}} H_{i}\left(x_{i}^{*}, u_{i}^{*}, p_{i+1}^{*}\right)=\nabla_{u_{i}}\left\{p_{i+1}^{*}{ }^{\prime} B_{i} u_{i}^{*}+\frac{1}{2} u_{i}^{* \prime} R_{i} u_{i}^{*}+\frac{1}{2} x_{i}^{* \prime} Q_{i} x_{i}^{*}\right\}=0
$$

or equivalently

$$
R_{i} u_{i}^{*}+B_{i}^{\prime} p_{i+1}^{*}=0, \quad i=0, \ldots, N-1
$$

from which

$$
u_{i}^{*}=-R_{i}^{-1} B_{i}^{\prime} p_{i+1}^{*}, \quad i=0, \ldots, N-1
$$

We will obtain a more convenient expression for $u_{i}^{*}$ by verifying a relation of the form

$$
p_{i+1}^{*}=K_{i+1} x_{i+1}^{*}
$$

where $K_{i+1}$ is a positive semidefinite symmetric matrix, which can be explicitly calculated. We show this relation by induction, first noting that it holds for $i=N-1$ with

$$
K_{N}=Q_{N}
$$

[cf. Eq. (9.42)]. Assume that it holds for some $i \leq N-1$. Then, by combining Eqs. (9.43) and (9.44), we have

$$
u_{i}^{*}=-R_{i}^{-1} B_{i}^{\prime} p_{i+1}^{*}=-R_{i}^{-1} B_{i}^{\prime} K_{i+1} x_{i+1}^{*}
$$

With the substitution $x_{i+1}^{*}=A_{i} x_{i}^{*}+B_{i} u_{i}^{*}$, this equation yields

$$
u_{i}^{*}=-R_{i}^{-1} B_{i}^{\prime} K_{i+1}\left(A_{i} x_{i}^{*}+B_{i} u_{i}^{*}\right)
$$

or equivalently, by solving for $u_{i}^{*}$,

$$
u_{i}^{*}=-\left(R_{i}+B_{i}^{\prime} K_{i+1} B_{i}\right)^{-1} B_{i}^{\prime} K_{i+1} A_{i} x_{i}^{*}, \quad i=0, \ldots, N-1
$$

We thus obtain

$$
x_{i+1}^{*}=A_{i} x_{i}^{*}+B_{i} u_{i}^{*}=A_{i} x_{i}^{*}-B_{i}\left(R_{i}+B_{i}^{\prime} K_{i+1} B_{i}\right)^{-1} B_{i}^{\prime} K_{i+1} A_{i} x_{i}^{*}
$$

Multiplying both sides with $A_{i}^{\prime} K_{i+1}$ and using Eq. (9.44), we see that

$$
A_{i}^{\prime} p_{i+1}^{*}=A_{i}^{\prime} K_{i+1} x_{i}^{*}-A_{i}^{\prime} K_{i+1} B_{i}\left(R_{i}+B_{i}^{\prime} K_{i+1}\right)^{-1} B_{i}^{\prime \prime} K_{i+1} B_{i} A_{i} x_{i}^{*}
$$

Adding $Q_{i} x_{i}^{*}$ to both sides, we obtain

$$
p_{i}^{*}=A_{i}^{\prime} p_{i+1}^{*}+Q_{i} x_{i}^{*}=K_{i} x_{i}^{*}
$$

where

$$
K_{i}=A_{i}^{\prime}\left(K_{i+1}-K_{i+1} B_{i}\left(R_{i}+B_{i}^{\prime} K_{i+1} B_{i}\right)^{-1} B_{i}^{\prime} K_{i+1}\right) A_{i}+Q_{i}
$$

thus completing the induction proof of Eq. (9.44).
Equation (9.47) is known as the Riccati equation. It generates the matrices $K_{1}, K_{2}, \ldots K_{N-1}$ from the terminal condition

$$
K_{N}=Q_{N}
$$

Given $K_{1}$, we may obtain $u_{0}^{*}$ using Eq. (9.46) and the given initial state $x_{0}$. Then, we may calculate $x_{1}^{*}$ from the system equation (9.39) and the vector $u_{1}^{*}$ from Eq. (9.46) using $K_{2}$. Similarly we may calculate $x_{1}^{*}, u_{2}^{*}, \ldots, u_{N-1}^{*}, x_{N}^{*}$.

---

# Gradient Methods for Optimal Control 

The application of gradient methods to unconstrained optimal control problems is straightforward in principle. For example the steepest descent method takes the form

$$
u_{i}^{k+1}=u_{i}^{k}-\alpha^{k} \nabla_{u_{i}} H_{i}\left(x_{i}^{k}, u_{i}^{k}, p_{i+1}^{k}\right), \quad i=0, \ldots, N-1
$$

where $H_{i}$ denotes the Hamiltonian function

$$
H_{i}\left(x_{i}, u_{i}, p_{i+1}\right)=g_{i}\left(x_{i}, u_{i}\right)+p_{i+1}^{\prime} f_{i}\left(x_{i}, u_{i}\right)
$$

$u^{k}=\left\{u_{0}^{k}, u_{1}^{k}, \ldots, u_{N-1}^{k}\right\}$ is the $k$ th control trajectory, $x^{k}=\left(x_{0}^{k}, x_{1}^{k}, \ldots, x_{N}^{k}\right)$ denotes the $k$ th state trajectory, and $p^{k}=\left(p_{1}{ }^{k} \ldots, p_{N}^{k}\right)$ denotes the $k$ th costate trajectory

$$
\begin{gathered}
p_{i}^{k}=\nabla_{x_{i}} H_{i}\left(x_{i}^{k}, u_{i}^{k}, p_{i+1}^{k}\right), \quad i=1, \ldots, N-1 \\
p_{N}^{k}=\nabla g_{N}\left(x_{N}^{k}\right)
\end{gathered}
$$

Thus, given $u^{k}$, one computes $x^{k}$ by forward propagation of the system equation, and then $p^{k}$ by backward propagation of the adjoint equation. Subsequently, the steepest descent iteration (9.48) is performed with the stepsize $\alpha^{k}$ chosen, for example, by the Armijo rule or some line minimization rule. An important point that guides the selection of the stepsize rule is that a gradient evaluation in optimal control problems with many control variables is usually not much more expensive than a function evaluation. Therefore, it is worth considering stepsize rules that aim at efficiency by using intelligently gradient calculations to save substantially on the number of function evaluations.

While steepest descent is simple, its convergence rate in optimal control problems is often very poor. This is particularly so when the underlying dynamic system tends to be unstable (see Exercise 9.1). For this reason scaling of the control variables is recommended if suitable scaling factors can be determined with reasonable effort. Otherwise, one should use either the conjugate gradient method or Newton's method.

## Conjugate Gradient Method for Optimal Control

The application of the conjugate gradient method to unconstrained optimal control problems is also straightforward. The search direction is a linear combination of the current gradient and the preceding search direction, as discussed in Section 1.6. The performance of the method can often be improved through the use of preconditioning. One possibility is to use a diagonal approximation to the Hessian matrix as preconditioning matrix, but in some cases, far more effective preconditioning schemes are possible: see Exercise 9.2 and [Ber74].

---

# Newton's Method for Optimal Control 

We first note that the next iterate of the pure form of Newton's method

$$
x^{k+1}=x^{k}-\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1} \nabla f\left(x^{k}\right)
$$

can be obtained by minimizing the second order Taylor series expansion of $f$ around $x^{k}$ given by

$$
f\left(x^{k}\right)+\nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)+\frac{1}{2}\left(x-x^{k}\right)^{\prime} \nabla^{2} f\left(x^{k}\right)\left(x-x^{k}\right)
$$

Carrying the argument one step further, assume that $\tilde{f}^{k}: \mathbb{R}^{n} \mapsto \mathbb{R}$ is a quadratic function such that

$$
\nabla \tilde{f}^{k}\left(x^{k}\right)=\nabla f\left(x^{k}\right), \quad \nabla^{2} \tilde{f}^{k}\left(x^{k}\right)=\nabla^{2} f\left(x^{k}\right)
$$

Then $\tilde{f}^{k}$ has the form

$$
\tilde{f}^{k}(x)=c+\nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)+\frac{1}{2}\left(x-x^{k}\right)^{\prime} \nabla^{2} f\left(x^{k}\right)\left(x-x^{k}\right)
$$

where $c$ is some constant, and hence the next iterate $x^{k+1}$ of Eq. (9.49) is obtained by minimizing the quadratic function $\tilde{f}^{k}$.

This viewpoint is particularly valuable in unconstrained optimal control problems, where it is impractical to evaluate explicitly the Hessian matrix $\nabla^{2} J(u)$ and find the control trajectory $u^{k+1}$ from

$$
u^{k+1}=u^{k}-\left(\nabla^{2} J\left(u^{k}\right)\right)^{-1} \nabla J\left(u^{k}\right)
$$

Instead, given $u^{k}$, we formulate a linear-quadratic optimal control problem with cost function $\tilde{J}_{k}(u)$ such that

$$
\nabla \tilde{J}_{k}\left(u^{k}\right)=\nabla J\left(u^{k}\right), \quad \nabla^{2} \tilde{J}_{k}\left(u^{k}\right)=\nabla^{2} J\left(u^{k}\right)
$$

The optimal solution $u^{k+1}$ of the linear quadratic problem can be obtained conveniently via a Riccati equation (cf. Example 8.3) and according to our previous discussion, it represents the next iterate of the pure form of Newton's method as applied to minimization of $J(u)$.

The linear quadratic problem that corresponds to a control trajectory $u^{k}=\left(u_{0}^{k}, u_{1}^{k}, \ldots u_{N-1}^{k}\right)$ and the associated states $x_{0}, x_{1}^{k}, \ldots x_{N}^{k}$, and costates $p_{1}^{k}, p_{2}^{k}, \ldots p_{N}^{k}$, involves a linearized version of the original system. It has the form

$$
\begin{aligned}
\operatorname{minimize} \tilde{J}_{k}(\delta u)= & \frac{1}{2} \delta x_{N}^{\prime} Q_{N} \delta x_{N}+a_{N}^{\prime} \delta x_{N}+\sum_{i=0}^{N-1} \frac{1}{2} \delta x_{i}^{\prime} Q_{i} \delta x_{i}+a_{i}^{\prime} \delta x_{i} \\
& +\sum_{i=0}^{N-1} \frac{1}{2} \delta u_{i}^{\prime} R_{i} \delta u_{i}+b_{i}^{\prime} \delta u_{i}+\sum_{i=1}^{N-1} \delta u_{i}^{\prime} M_{i} \delta x_{i}
\end{aligned}
$$

subject to $\delta x_{0}=0, \quad \delta x_{i+1}=A_{i} \delta x_{i}+B_{i} \delta u_{i}, \quad i=1, \ldots, N-1$,

---

where

$$
\begin{gathered}
A_{i}=\nabla_{x_{i}} f_{i}^{\prime}, \quad B_{i}=\nabla_{u_{i}} f_{i}^{\prime}, \quad i=0, \ldots, N-1 \\
Q_{N}=\nabla^{2} g_{N}, \quad a_{N}=\nabla g_{N} \\
Q_{i}=\nabla_{x_{i} x_{i}}^{2} H_{i}, \quad a_{i}=\nabla_{x_{i}} g_{i}, \quad i=1, \ldots, N-1 \\
R_{i}=\nabla_{u_{i} u_{i}}^{2} H_{i}, \quad b_{i}=\nabla_{u_{i}} g_{i}, \quad i=0, \ldots, N-1 \\
M_{i}=\nabla_{u_{i} x_{i}}^{2} H_{i}, \quad i=0, \ldots, N-1 \\
H_{i}\left(x_{i}, u_{i}, p_{i+1}^{k}\right)=g_{i}\left(x_{i}, u_{i}\right)+p_{i+1}^{k}{ }^{\prime} f_{i}\left(x_{i}, u_{i}\right)
\end{gathered}
$$

The partial derivatives appearing above are evaluated along the current control and state trajectories

$$
u^{k}=\left(u_{0}^{k}, u_{1}^{k}, \ldots, u_{N-1}^{k}\right), \quad x^{k}=\left(x_{0}, x_{1}^{k}, \ldots, x_{N}^{k}\right)
$$

and the corresponding costate vectors $p_{1}^{k}, p_{x}^{k}, \ldots, p_{N}^{k}$. The problem (9.51) involves the variations

$$
\begin{aligned}
& \delta x_{i}=x_{i}-\tilde{x}_{i}^{k} \\
& \delta u_{i}=u_{i}-u_{i}^{k}
\end{aligned}
$$

from the nominal trajectories $u^{k}, \tilde{x}^{k}$, where $\tilde{x}^{k}$ is the trajectory corresponding to the linearized system

$$
\begin{gathered}
\tilde{x}_{i+1}^{k}=A_{i} \tilde{x}_{i}^{k}+B_{i} u_{i}^{k}, \quad i=0, \ldots, N-1 \\
\tilde{x}_{0}^{k}=x_{0}: \text { given initial state. }
\end{gathered}
$$

It is straightforward to verify that the cost function $J(u)$ for our original problem has a gradient $\nabla J\left(u^{k}\right)$ and Hessian matrix $\nabla^{2} J\left(u^{k}\right)$ evaluated at $u^{k}$, which are equal to the corresponding gradient and Hessian matrix $\nabla \hat{J}_{k}\left(\delta u^{k}\right), \nabla^{2} \hat{J}_{k}\left(\delta u^{k}\right)$ of the problem (9.51)-(9.52). This can be shown by comparing the expressions for the gradient and Hessian of the two problems given earlier [cf. Eqs. (9.29)-(9.31)]. Hence, according to our earlier discussion, the solution of problem (9.51)-(9.52) yields the next iterate of Newton's method.

We now turn to the linear-quadratic problem (9.51)-(9.52). It is possible to show using the necessary conditions for optimality, similarly as in Example 8.3, that the solution of this problem (provided it exists and is unique) may be obtained in feedback form as follows:

$$
\begin{gathered}
\delta u_{i}^{k}=-\left(R_{i}+B_{i}^{\prime} K_{i+1} B_{i}\right)^{-1}\left(\left(M_{i}+B_{i}^{\prime} K_{i+1} A_{i}\right) \delta x_{i}^{k}+b_{i}+B_{i}^{\prime} \lambda_{i+1}\right) \\
i=0, \ldots, N-1
\end{gathered}
$$

where

$$
\delta x_{0}^{k}=0
$$

---

$$
\delta x_{i+1}^{k}=A_{i} \delta x_{i}^{k}+B_{i} \delta u_{i}^{k}, \quad i=0, \ldots, N-1
$$

the matrices $K_{1}, K_{2}, \ldots, K_{N}$ are given recursively by the equation

$$
\begin{gathered}
K_{N}=Q_{N} \\
K_{i}=A_{i}^{\prime} K_{i+1} A_{i}+Q_{i}-\left(B_{i}^{\prime} K_{i+1} A_{i}+M_{i}\right)^{\prime}\left(R_{i}+B_{i}^{\prime} K_{i+1} B_{i}\right)^{-1} \\
\left(B_{i}^{\prime} K_{i+1} A_{i}+M_{i}\right), \quad i=1, \ldots, N-1
\end{gathered}
$$

and the vectors $\lambda_{1}, \lambda_{2}, \ldots, \lambda_{N}$ are given recursively by

$$
\begin{gathered}
\lambda_{N}=a_{N} \\
\lambda_{i}=a_{i}+A_{i}^{\prime} \lambda_{i+1}-A_{i}^{\prime} K_{i+1} B_{i}\left(R_{i}+B_{i}^{\prime} K_{i+1} B_{i}\right)^{-1}\left(b_{i}+B_{i}^{\prime} \lambda_{i+1}\right)
\end{gathered}
$$

The next iterate of the pure form of Newton's method is obtained from Eq. (9.55) and if a stepsize $\alpha^{k}$ is introduced, then we obtain the iteration

$$
u_{i}^{k+1}=u_{i}^{k}+\alpha^{k} \delta u_{i}^{k}, \quad i=0, \ldots, N-1
$$

with

$$
\delta u_{i}^{k}=-\left(R_{i}+B_{i}^{\prime} K_{i+1} B_{i}\right)^{-1}\left(\left(M_{i}+B_{i}^{\prime} K_{i+1} A_{i}\right) \delta x_{i}^{k}+b_{i}+B_{i}^{\prime} \lambda_{i+1}\right)
$$

The computations of the typical iteration of Newton's method are carried out in the following sequence.

# kth Iteration of Newton's Method 

Step 1: Given the current control trajectory $u^{k}$ the corresponding state trajectory $x^{k}$ is calculated.
Step 2: The solution $K_{1}, \ldots, K_{N}$ of Eq. (9.57) and the solution $\lambda_{1}, \ldots, \lambda_{N}$ of Eq. (9.58) are calculated (backwards, from the terminal conditions) together with all the necessary data for calculation of the optimal control in equation (9.55). This is done by calculating simultaneously the costate vectors $p_{1}^{k}, \ldots, p_{N}^{k}$, which are needed to evaluate the second derivatives of the Hamiltonian (i.e., $Q_{i}, R_{i}, M_{i}$ ).

Step 3: The Newton direction $\delta u^{k}$ is generated as follows:

1) Calculate $\delta u_{0}^{k}$ using Eq. (9.60),

$$
\delta u_{0}^{k}=-\left(R_{0}+B_{0}^{\prime} K_{1} B_{0}\right)^{-1}\left(b_{0}+B_{0}^{\prime} \lambda_{1}\right)
$$

2) Calculate $\delta x_{1}^{k}$ using Eq. (9.56) and the fact $\delta x_{0}^{k+1}=0$.
3) Calculate $\delta u_{1}^{k}$ using $\delta x_{1}^{k}$ and Eq. (9.60).
4) Continue similarly to compute $\delta u_{i}^{k}$ and $\delta x_{i}^{k}$ for all $i$.

---

Step 4: A stepsize $\alpha^{k}$ is obtained using some stepsize rule; the next control trajectory is

$$
u_{i}^{k+1}=u_{i}^{k}+\alpha^{k} \delta u_{i}^{k}, \quad i=0, \ldots, N-1
$$

In the computation of the solution of the Riecati equation it is necessary that $\left(R_{i}+B_{i}^{\prime} K_{i+1} B_{i}\right)^{-1}$ exist and be positive definite. If any of these inverses is not positive definite, it can be seen that $\nabla^{2} J\left(u^{k}\right)$ is not positive definite and hence the Newton iteration must be modified. Under these circumstances one should replace $R_{i}$ by $\left(R_{i}+E_{i}\right)$ where $E_{i}$ is a positive definite diagonal matrix such that $\left(R_{i}+E_{i}+B_{i}^{\prime} K_{i+1} B_{i}\right)$ is positive definite, as discussed in Section 1.4. In a neighborhood of a local minimum $u^{*}$ with $\nabla^{2} J\left(u^{*}\right)>0$, the solution of the Riecati equation exists and the Newton iteration can be carried out as described.

It is interesting to note that our analysis yields as a byproduct a sufficient condition for optimality of a control trajectory $u^{*}$. Clearly $u^{*}$ will be optimal if the corresponding linear-quadratic problem (9.51) with all derivatives evaluated at $u^{*}$ has $\delta u=0$ as an optimal solution. This can be guaranteed when

$$
\begin{gathered}
\nabla_{u_{i}} H_{i}=0, \quad i=0, \ldots, N-1, \quad \text { (1st order condition), } \\
R^{*}+B_{i}^{* \prime} K_{i+1}^{*} B_{i}^{*}>0, \quad i=0, \ldots, N-1, \quad \text { (2nd order condition) }
\end{gathered}
$$

where

$$
R_{i}^{*}=\nabla_{u_{i}=}^{2} H_{i}, \quad B_{i}^{*}=\nabla_{u_{i}} f_{i}^{\prime}, \quad i=0, \ldots, N-1
$$

$K_{i+1}^{*}$ is given by the Riecati equation (9.57) corresponding to $u^{*}$, and all derivatives above are evaluated along $u^{*}$ and the corresponding state and costate trajectories.

Finally, we mention a variation of Newton's method that offers some computational advantages. It uses in place of Eqs. (9.59) and (9.60) the following two equations:

$$
\begin{aligned}
u_{i}^{k+1}=u_{i}^{k}-\alpha^{k}\left(R_{i}+B_{i}^{\prime} K_{i+1} B_{i}\right)^{-1} & \left(\left(M_{i}+B_{i}^{\prime} K_{i-1} A_{i}\right)\left(x_{i}^{k+1}-x_{i}^{k}\right)\right. \\
& \left.+A_{i}+B_{i}^{\prime} \lambda_{i+1}\right)
\end{aligned}
$$

where

$$
x_{i}^{k+1}-x_{i}^{k}=f_{i-1}\left(x_{i-1}^{k+1} \cdot u_{i-1}^{k+1}\right)-f_{i-1}\left(x_{i-1}^{k}, u_{i-1}^{k}\right)
$$

In other words, the first order variation $\delta x_{i}^{k}$, which is obtained from the linearized system equation (9.56), is replaced in Eq. (9.60) by the actual variation $\left(x_{i}^{k+1}-x_{i}^{k}\right)$ of Eq. (9.62). Since $\left(x_{i}^{k+1}-x_{i}^{k}\right)$ and $\delta x_{i}^{k}$ differ by terms which are of second or higher order, the two methods are asymptotically equivalent. The advantage of using Eqs. (9.61) and (9.62) is that the matrices $A_{i}, B_{i}$ need not be stored (or recalculated) for use in Eq. (9.56). In addition, the trajectory $x^{k+1}$ corresponding to $u^{k+1}$ may be used in the next iteration, if the pure form of Newton's method is used or the stepsize $\alpha^{k}$ turns out to be equal to one.

---

# E X E R C I S E S 

## 9.1

Calculate the condition number of the optimal control problem involving the quadratic cost function $x_{N}^{2}+\sum_{i=0}^{N-1} u_{i}^{2}$ and the scalar system $x_{i+1}=a x_{i}+u_{i}$, where $a>1$. Show that it tends to $\infty$ as $N \rightarrow \infty$.

## $9.2[\operatorname{Ber} 74]$

Consider the special case of the linear quadratic problem of Eqs. (9.39) and (9.40) where $Q_{i}=0$ for all $i=0, \ldots, N-1$. Derive a preconditioned conjugate gradient method that converges in at most $n+1$ steps. Hint: Consider the structure of the Hessian matrix and use the result of Exercise 6.2 of Section 1.6 .

### 1.10 SOME PRACTICAL GUIDELINES

Practical nonlinear programming problems can be very challenging. They require an iterative process that can be slow and may in the end lead to just a local minimum. It is thus useful to follow a few basic practical guidelines. Keep in mind, however, that in nonlinear programming, just as there is no foolproof method, there is also no foolproof advice.

## Problem Formulation

In the real world, optimization problems seldom come neatly formulated as mathematical problems. Thus, usually the first step is to select the cost function and to delineate the constraints. One potential difficulty here is that there may be multiple and possibly competing objectives; for example in an engineering design problem one may simultaneously want to minimize cost and maximize efficiency. In this case one usually encodes all objectives except one in the constraints and embodies the remaining objective in the cost function. However, the division between objective function and constraints may be arbitrary and one may have to reconsider the formulation after evaluating the results of the optimization. It is important here to properly correlate the problem formulation with the objectives of the investigation. For example if only an approximate answer is required, there is no point in constructing a highly detailed model.

---

The most critical task in problem formulation is to capture the realism of the practical context while obtaining an analytically or computationally tractable model. This is where experience and insight into both the application and the methodology is very important. In particular, it is important to be able to recognize the models that are easily solvable and those that are essentially impossible to solve. Two rules of thumb here are:
(a) Convex cost and constraint sets are preferable to nonconvex ones.
(b) Linear programming problems are easier than nonlinear programming problems, which are in turn easier than integer/discrete combinatorial problems.
While (a) is generally true because of the lack of local minima and other spurious stationary points in convex problems, there are many exceptions to (b). For example there is a widespread practice of replacing nonlinear functions by piecewise linear ones that lead to linear programming formulations, so that linear programming codes can be used. Unfortunately this may lead to vastly increased dimensionality and loss of insight (solutions to continuous problems often have elegant features that are lost in discrete approximations). Another example along similar lines arises when essentially combinatorial problems are formulated as highly nonconvex nonlinear programming problems whose local minima correspond to the feasible solutions of the original problem. For a somewhat absurd example, think of discarding an integrality constraint on a variable $x_{i}$ and adding to the cost function the (differentiable) term $c\left(x_{i}-n\left(x_{i}\right)\right)^{2}$ where $n\left(x_{i}\right)$ is the nearest integer to $x_{i}$ and $c$ is a very large penalty parameter. Generally, when a convex problem formulation seems impossible and the presence of many local minima is an important concern, one may wish to bring to bear the methods of global optimization, for which we refer to the specialized literature (see e.g., [FlP95], [Flo95], [HPT95], [PaR87]).

We note, however, that depending on the practical situation, there is often merit in replacing integer variables with continuous variables and using a nonlinear programming formulation to obtain a noninteger solution that can then be rounded to integer. For example, consider a variable such as number of telephone lines to establish between two points in a communication network. This number is naturally discrete but if its range is in the hundreds, little will be lost by replacing it by a continuous variable. On the other hand, other variables may naturally be discrete-valued, such as a $\{0,1\}$-valued variable modeling whether route A or route B is used to send a message in a communication network. A continuous approximation of such a variable may be meaningless.

Another important consideration is whether the cost function is or is not differentiable. The most powerful methods in nonlinear programming (e.g. the ones discussed in this chapter) require a once or twice differentiable cost function. This motivates the smoothing of whatever natural

---

nondifferentiabilities may exist in the cost. An important special case is when nondifferentiable terms of the form $\max \left\{f_{1}(x), f_{2}(x)\right\}$ are replaced in the cost function and/or the constraints by a smooth approximation (see Fig. 1.10.1 and Exercise 4.6 in Section 6.4). On the other hand there are situations, particularly arising in the context of duality where nondifferentiable cost functions are unavoidable and the use of a nondifferentiable optimization method is naturally suited to the character of the problem (see Chapter 6). Note that the implementation of gradient methods can be enhanced with the use of automatic differentiation programs that compute first and second derivatives from user-supplied programs that compute only function values; see [Gri89], [Gil92].

# Scaling 

When selecting the optimization variables in a given problem it is important to pay attention to their natural relative magnitudes, and scale them so that their values are neither too large nor too small. This is useful to control roundoff error and it is also often helpful for improving the condition number of the problem, and the natural convergence rate of steepest descent and other algorithms. However, one should also pay attention to the range of the variables as well, and use a translation to somewhere in the middle of that range. In particular, suppose that the optimization variable $x_{i}$ is expected to take values in the range $\left[\alpha_{i}, \beta_{i}\right]$. Then it is typically helpful to use the transformation of variables

$$
y_{i}=\frac{x_{i}-\left(\alpha_{i}+\beta_{i}\right) / 2}{\left(\beta_{i}-\alpha_{i}\right) / 2}
$$

so that the new variable $y_{i}$ takes values in the range $[-1,1]$. Note that this preliminary scaling of the variables complements but is not a substitute for the iteration-dependent scaling that we discussed in Section 1.3 and which is inherent in Newton-like methods.

Some least squares problems are special because the variables can be scaled automatically by scaling some of the coefficients in the data blocks. As an example consider the linear least squares problem

$$
\begin{aligned}
& \operatorname{minimize} \sum_{i=1}^{m}\left(y_{i}-a_{i}^{\prime} x\right)^{2} \\
& \text { subject to } x \in \Re^{n}
\end{aligned}
$$

where $y_{i}$ are given scalars and $a_{i}$ are given vectors in $\Re^{n}$ with coordinates $a_{i j}, j=1, \ldots, n$. The $j$ th diagonal element of the Hessian matrix is

$$
d_{j}=\sum_{i=1}^{m} a_{i j}^{2}
$$

---

![[chapter_21_p168_img31.jpeg]]

Figure 1.10.1. Smoothing of nondifferentiable terms. The function $f: \Re^{n} \rightarrow \Re$ given by

$$
f(x)=\max \left\{f_{1}(x), f_{2}(x)\right\}=f_{1}(x)+\max \left\{0, f_{2}(x)-f_{1}(x)\right\}
$$

is replaced by

$$
f_{1}(x)+\tilde{\gamma}\left(f_{2}(x)-f_{1}(x), \lambda, c\right)
$$

where $\tilde{\gamma}(t, \lambda, c)$ is a differentiable function that approximates $\gamma(t)=\max \{0, t\}$. This function depends on a parameter $\lambda \in[0,1]$ and a parameter $c>0$. It has the form

$$
\tilde{\gamma}(t, \lambda, c)= \begin{cases}t-\frac{(1-\lambda)^{2}}{2 c} & \text { if } \frac{1-\lambda}{c} \leq t \\ \lambda t+\frac{c}{2} t^{2} & \text { if }-\frac{\lambda}{c} \leq t \leq \frac{1-\lambda}{c} \\ -\frac{\lambda^{2}}{2 c} & \text { if } t \leq-\frac{\lambda}{c}\end{cases}
$$

shown in the figure. The accuracy of the approximation increases as $c$ increases. The scalar $\lambda$ plays the role of a Lagrange multiplier; see Exercise 4.6 in Section 5.4. This smoothing method can be generalized for the case

$$
f(x)=\max \left\{f_{1}(x), \ldots, f_{m}(x)\right\}
$$

see [Ber75c], [Ber77], [Geo77], [Ber82a], [Pap81], [Pol79] for further analysis. An alternative smooth approximation of $f$ is the function

$$
\frac{1}{c} \ln \left\{\sum_{i=1}^{m} \lambda_{i} e^{c f_{i}(x)}\right\}
$$

where $c>0$ and $\lambda_{i}$ are positive numbers with $\sum_{i=1}^{m} \lambda_{i}=1$ (see [Ber82a], Section 5.1.3).

Thus reasonable scaling of the variables is obtained by multiplying for each $j$, all coefficients $a_{i j}$ with a common scalar $s_{j}$ so that they all lie in the range $[-1,1]$.

---

# Method Selection 

An important consideration here is whether a one-time solution of a given problem is sought or whether the problem (with data variations) will be solved many times on a production basis. In the former case any method that will do the job is sufficient, so one should aim for speed of code development, perhaps through the use of an existing optimization package. In the latter case, efficiency and accuracy of solution become important, so careful experimentation with a variety of methods may be appropriate in order to obtain a method that solves the problem fast and exploits its special structure.

Throughout this chapter we emphasized the fast convergence rate of Newton-like methods when derivatives are available and the problem is nonsingular. It is worth repeating, however, that for singular and difficult problems, simpler steepest descent-like methods with simple scaling may be more effective. Generally, in lack of a clear choice it is best to start with a simple method that can be easily coded and understood, and be prepared to proceed with more sophisticated choices.

Once a type of method is selected, it must be tuned to the problem at hand. In particular, stepsize rule parameters and termination criteria must be chosen. This often requires trial and error, particularly for first order methods, for a diminishing or constant stepsize, and for difficult (e.g., singular) problems. Termination criteria should be scale-free to the extent possible. For example it is preferable to use as termination tests $\left\|x^{k}-x^{k-1}\right\| \leq \epsilon\left\|x^{k}\right\|$ or $\left\|\nabla f\left(x^{k}\right)\right\| \leq \epsilon\left\|\nabla f\left(x^{0}\right)\right\|$ rather than $\left\|x^{k}-x^{k-1}\right\| \leq \epsilon$ or $\left\|\nabla f\left(x^{k}\right)\right\| \leq \epsilon$, respectively.

## Validation

Finally one must strive to be convinced that the answers obtained from the computation are reasonable approximations to a (global) minimum. There are no systematic methods for this, but a few heuristics tailored to the problem at hand are adequate in most cases. In particular, suppose that some special case or variation of the given problem has a known optimal solution. Then it is reassuring if the coded method succeeds in finding this solution. Similarly, if a lower bound to the optimal value is known, it is reassuring when the algorithm comes close to achieving this lower bound. For example, in least squares problems, a small cost function value is a strong indication of success.

Another frequently used technique is to start the algorithm from widely varying starting points and see if it will give a better answer. This is particularly recommended in the case of multiple local minima. Finally, one may vary the parameters of the algorithm and also try alternative algorithms to see if substantial improvements can be obtained.

---

# 1.11 NOTES AND SOURCES 

Section 1.2: The steepest descent method dates to Cauchy [Cau47], who attributes to Newton the unity stepsize version of what we call Newton's method. The modern theory of gradient methods evolved in the 60's, when on one hand, several practical stepsize rules were proposed starting with the work of Goldstein [Gol62], [Gol64], and on the other hand, general methods of convergence analysis were formulated in the works of Zangwill [Zan69], Polak [Pol71], Daniel [Dan71], and Ortega and Rheinboldt [OrR70]. The Armijo stepsize rule [Arm66] is the most popular of a broad variety of rules that enforce descent without requiring a full line minimization. Stepsize rules that do not enforce strict descent at each iteration, while aiming at faster convergence, are discussed in [GLL91], and also in [BaB88] and [Ray93]. The capture theorem (Prop. 1.2.5) was formulated and proved by the author for the case of a nonsingular local minimum in [Ber82a], and was extended to the form given here by Dunn [Dun93c]. Gradient methods with errors and diminishing stepsize are discussed in [Pol87], which summarizes much earlier work on the subject. Parallel and asynchronous versions of these methods converge under very weak conditions [TBA86], [BeT89].

Section 1.3: For further discussion on various measures of rate of convergence, see [OrR70], [Ber82a], and [BaD93]. The convergence rate of steepest descent with line minimization was analyzed by Kantorovich [Kan45]. The case of a constant stepsize was analyzed in [Gol64] and [LeP65]. For analysis of the convergence rate of steepest descent for singular problems see [Dun81], [Dun87], and [Pol87].

Section 1.4: The local convergence of Newton's method was analyzed by Kantorovich [Kan49]. For an analysis of the case where the method converges to a singular point, see [DeK80], [DKK83], and [HuD84]. The modification to extend the region of convergence of Newton's method by modifying the Cholesky factorization is described in [GiM74]: see also [GMW81]. The use of a trust region has been discussed by [MoS83]. Extensive accounts of various aspects of Newton-like methods are given in [Gol67], [GMW81], [DeS83], and [Naz91]. The truncated Newton method is discussed in [DES82], [Nas85], and [NaS89].

Section 1.5: Incremental gradient methods for linear least squares problems are attributed to [WiH59]. For recent analyses, see [Gai94], [Gri94], [LuT94a], [MaS94]. Proposition 1.5.1 is due to [Luo91] and stems from an earlier result of [Koh74]. The incremental Gauss-Newton method was proposed in [Dav76] on an empirical basis. The development given here in terms of the Extended Kalman Filter and the corresponding convergence analysis are new [Ber94]. In the case where $\lambda<1$, and for some $x^{*}$ we have $g\left(x^{*}\right)=0$ and $\nabla g\left(x^{*}\right)$ : nonsingular, the method has been shown to converge to $x^{*}$ linearly with convergence ratio $\lambda$ when started sufficiently

---

close to $x^{*}$ [Pap82]. For a description of methods for neural network training problems, see [Hay94]. The ill-conditioned character of these problems is discussed in [SBC93].

Section 1.6: Conjugate direction methods are due to Hestenes and Stiefel [HeS52]. The convergence rate of the conjugate gradient method is discussed in [FaF63], and also in [Lue84], which we follow in our discussion. The conjugate gradient method with the formula (6.36) is known as the Poljak-Polak-Ribiere method and its convergence was analyzed in [Pol69a] and [PoR69]. An extensive discussion of conjugate direction methods can be found in the monograph [Hes80].
Section 1.7: The DFP method is due to Davidon [Dav59]. It was popularized primarily through the work of Fletcher and Powell [FlP63]. The BFGS method is independently due to [Bro70], [Fle70a], [Gol70], [Sha70]. Surveys and extensive discussions of Quasi-Newton methods include [DeM77], [DeS83], and [Naz91].

Section 1.8: For further material on direct search methods, see [GiM74] and $[\mathrm{Avr} 76]$.

Section 1.9: An extensive reference on optimality conditions for discretetime optimal control is [CCP70]; for computational methods see [Pol71] and also the survey [Pol73]. The implementation of Newton's method given here is due to [DuB89] and was also independently derived in [PaM89] (see [Mit66] for a continuous-time version of this algorithm). An alternative approach to discrete-time optimal control is based on dynamic programming; see e.g. [Ber95a]. Discrete-time optimal control problems often arise from discretization of their continuous-time counterparts. The associated issues are discussed in [Cul71].

Section 1.10: A great deal of material on the implementation of nonlinear programming methods is given in [GMW81]. A description of available software packages is provided in [MoW93].

---

.