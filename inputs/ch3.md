#  Lagrange Multiplier Theory 

## Contents

3.1. Necessary Conditions for Equality Constraints ..... p. 254
3.2. Sufficient Conditions and Sensitivity Analysis ..... p. 271
3.3. Inequality Constraints ..... p. 282
3.4. Linear Constraints* ..... p. 292
3.5. Notes and Sources ..... p. 310

---

The methods I set forth require neither constructions nor geometric or mechanical considerations.
They require only algebraic operations subject to a systematic and uniform course.

# Lagrange 

The constraint set of an optimization problem is usually specified in terms of equality and inequality constraints. If we take into account this structure, we obtain a sophisticated collection of optimality conditions, involving some auxiliary variables called Lagrange multipliers. These variables facilitate the characterization of optimal solutions, but also provide valuable sensitivity information, quantifying to first order the variation of the optimal cost caused by variations in problem data.

Lagrange multipliers will be central in almost all of the subsequent material in the book. In Chapter 4, we will see how they can be used within computational methods. In Chapter 5, we will see how they can be viewed as the optimization variables of auxiliary optimization problems, called dual problems. The computational solution of dual problems will be the subject of Chapter 6. The present chapter also contains an introductory discussion of duality.

### 3.1 NECESSARY CONDITIONS FOR EQUALITY CONSTRAINTS

In this section we consider problems with equality constraints of the form

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } h_{i}(x)=0, \quad i=1, \ldots, m
\end{aligned}
$$

We assume that $f: \Re^{n} \mapsto \Re, h_{i}: \Re^{n} \mapsto \Re, i=1, \ldots, m$, are continuously differentiable functions. All the necessary and the sufficient conditions of this chapter relating to a local minimum can also be shown to hold if $f$ and $h_{i}$ are defined and are continuously differentiable in just an open set containing the local minimum. The proofs are essentially identical to those given here.

For notational convenience, we introduce the constraint function $h$ : $\Re^{n} \mapsto \Re^{m}$, where

$$
h=\left(h_{1}, \ldots, h_{m}\right)
$$

We can then write the constraints in the more compact form

$$
h(x)=0
$$

---

Our basic Lagrange multiplier theorem states that for a given local minimum $x^{*}$, there exist scalars $\lambda_{1}, \ldots, \lambda_{m}$, called Lagrange multipliers, such that

$$
\nabla f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i} \nabla h_{i}\left(x^{*}\right)=0
$$

There are two ways to interpret this equation:
(a) The cost gradient $\nabla f\left(x^{*}\right)$ belongs to the subspace spanned by the constraint gradients. The example of Fig. 3.1.1 illustrates this interpretation.
(b) The cost gradient $\nabla f\left(x^{*}\right)$ is orthogonal to the subspace of first order feasible variations

$$
V\left(x^{*}\right)=\left\{\Delta x \mid \nabla h_{i}\left(x^{*}\right)^{\prime} \Delta x=0, i=1, \ldots, m\right\}
$$

This is the subspace of variations $\Delta x$ for which the vector $x=x^{*}+\Delta x$ satisfies the constraint $h(x)=0$ up to first order. Thus, according to the Lagrange multiplier condition of Eq. (1.2), at the local minimum $x^{*}$, the first order cost variation $\nabla f\left(x^{*}\right)^{\prime} \Delta x$ is zero for all variations $\Delta x$ in this subspace. This is an analogous statement to the "zero gradient condition" $\nabla f\left(x^{*}\right)=0$ of unconstrained optimization.

The formal statement of the main Lagrange multiplier theorem is as follows:

Proposition 3.1.1: (Lagrange Multiplier Theorem - Necessary Conditions) Let $x^{*}$ be a local minimum of $f$ subject to $h(x)=$ 0 , and assume that the constraint gradients $\nabla h_{1}\left(x^{*}\right), \ldots, \nabla h_{m}\left(x^{*}\right)$ are linearly independent. Then there exists a unique vector $\lambda^{*}=$ $\left(\lambda_{1}^{*}, \ldots, \lambda_{m}^{*}\right)$, called a Lagrange multiplier vector, such that

$$
\nabla f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla h_{i}\left(x^{*}\right)=0
\tag{1.3}
$$

^8eef07

If in addition $f$ and $h$ are twice continuously differentiable, we have

$$
y^{\prime}\left(\nabla^{2} f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla^{2} h_{i}\left(x^{*}\right)\right) y \geq 0, \quad \text { for all } y \in V\left(x^{*}\right)
\tag{1.4}
$$

^2a7ec6

where $V\left(x^{*}\right)$ is the subspace of first order feasible variations

$$
V\left(x^{*}\right)=\left\{y \mid \nabla h_{i}\left(x^{*}\right)^{\prime} y=0, i=1, \ldots, m\right\}
\tag{1.5}
$$

---

![[ch3_and_4_p4_img1.jpeg]]

Figure 3.1.1. Illustration of the Lagrange multiplier condition (1.1) for the problem

$$
\begin{aligned}
& \operatorname{minimize} x_{1}+x_{2} \\
& \text { subject to } x_{1}^{2}+x_{2}^{2}=2
\end{aligned}
$$

At the local minimum $x^{*}=(-1,-1)$, the cost gradient $\nabla f\left(x^{*}\right)$ is normal to the constraint surface and is therefore, collinear with the constraint gradient $\nabla h\left(x^{*}\right)=(-2,-2)$. The Lagrange multiplier is $\lambda=1 / 2$.

For easy reference, a feasible vector $x$ for which the constraint gradients $\nabla h_{1}(x), \ldots, \nabla h_{m}(x)$ are linearly independent will be called **regular**. We will see later that *there may not exist Lagrange multipliers for a local minimum, which is not regular*.

We will provide two different proofs of the Lagrange multiplier theorem, each providing important insights. These proofs are based on transforming the constrained problem to an unconstrained one, but in different ways. The constrained first and second order necessary conditions are obtained by applying the corresponding unconstrained conditions to the appropriate unconstrained problem. The approaches are:
(a) **The elimination approach**. Here we view the constraints as a system of $m$ equations with $n$ unknowns, and we express $m$ of the variables in terms of the remaining $n-m$, thereby reducing the problem to an unconstrained problem. By then applying the corresponding first and second order necessary conditions for unconstrained minima, the

---

Lagrange multiplier theorem follows. This approach requires the use of the **implicit function theorem** ([[Mathematics/Optimization/Nonlinear Optimization/Textbook/Appendix#^aef75a|Prop. A. 25]] in Appendix A) but is otherwise simple and insightful. Its extension, however, to inequality constraints is cumbersome.
(b) **The penalty approach**. Here we disregard the constraints, while adding to the cost a high penalty for violating them. By writing the necessary conditions for the "penalized" unconstrained problems, and by passing to the limit as the penalty increases, we obtain the Lagrange multiplier theorem. This approach does not require the use of the implicit function theorem. It is simple, and applies to inequality constraints, as will be seen in Section 3.3.

# 3.1.1 The Elimination Approach 

We introduce the elimination approach by considering first the easier case where the constraints are linear.

## Example 1.1 (Lagrange Multipliers for Linear Constraints)

Consider the problem

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } A x=b
\end{aligned}
\tag{1.6}
$$

^ff86fd

where $A$ is an $m \times n$ matrix with linearly independent rows and $b \in \Re^{m}$ is a given vector. By rearranging the coordinates of $x$ if necessary, [[Note 3 - Partitioning matrix A|we may assume that]] the first $m$ columns of $A$ are linearly independent, so that $A$ can be partitioned as

$$
A=\left(\begin{array}{ll}
B & R
\end{array}\right)
$$

where $B$ is an invertible $m \times m$ matrix and $R$ is an $m \times(n-m)$ matrix. We also partition $x$ as

$$
x=\binom{x_{B}}{x_{R}}
$$

where $x_{B} \in \Re^{m}$ and $x_{R} \in \Re^{n-m}$. 

>*Austin Note*: We want to partition this way so that we can eliminate a variable. This will allow us to isolate $x_B$ by taking inverse of $B$

We can then write problem (1.6) as

$$
\begin{aligned}
& \operatorname{minimize} f\left(x_{B}, x_{R}\right) \\
& \text { subject to } B x_{B}+R x_{R}=b
\end{aligned}
\tag{1.7}
$$

By using the constraint equation to express $x_{B}$ in terms of $x_{R}$ as

$$
x_{B}=B^{-1}\left(b-R x_{R}\right)
$$

and by substitution in the cost function, we can convert problem (1.7) into the unconstrained optimization problem

$$
\begin{aligned}
& \operatorname{minimize} F\left(x_{R}\right) \equiv f\left(B^{-1}\left(b-R x_{R}\right), x_{R}\right) \\
& \text { subject to } x_{R} \in \Re^{n-m}
\end{aligned}
\tag{1.8}
$$

^1b0173

---

If $\left(x_{B}^{*}, x_{R}^{*}\right)$ is a local minimum of the constrained problem (1.6), then $x_{R}^{*}$ is [[Note 4 - Why optimality transfers|an unconstrained local minimum of the "reduced" cost function]] $F$, so we have

$$
0=\nabla F\left(x_{R}^{*}\right)=-R^{\prime}\left(B^{\prime}\right)^{-1} \nabla_{B} f\left(x^{*}\right)+\nabla_{R} f\left(x^{*}\right)
\tag{1.9}
$$

^f44962

>[[Note 5 - Derivation of grad F for reduced cost function]]

where $\nabla_{B} f$ and $\nabla_{R} f$ denote the gradients of $f$ with respect to $x_{B}$ and $x_{R}$, respectively. By defining

$$
\lambda^{*}=-\left(B^{\prime}\right)^{-1} \nabla_{B} f\left(x^{*}\right)
\tag{1.10}
$$

Eq. (1.9) is written as

$$
\nabla_{R} f\left(x^{*}\right)+R^{\prime} \lambda^{*}=0
\tag{1.11}
$$

while Eq. (1.10) is written as

$$
\nabla_{B} f\left(x^{*}\right)+B^{\prime} \lambda^{*}=0
\tag{1.12}
$$

Equations (1.11) and (1.12) can be combined as

$$
\nabla f\left(x^{*}\right)+A^{\prime} \lambda^{*}=0
\tag{1.13}
$$

>[[Note 6 - Lagrange multiplier combination step]]

which is the Lagrange multiplier condition of [[Chapter 3#^8eef07|Eq. (1.3)]] specialized to the case of the linearly constrained problem [[Chapter 3#^ff86fd|(1.6)]]. Note that the vector $\lambda^{*}$ satisfying this condition is unique because of the linear independence of the columns of $A^{\prime}$.

We next show the second order necessary condition [[Chapter 3#^2a7ec6|(1.4)]], by showing that it is equivalent to the second order unconstrained necessary condition

$$
0 \leq d^{\prime} \nabla^{2} F\left(x_{R}^{*}\right) d, \quad \forall d \in \Re^{n-m}
\tag{1.14}
$$

We have using [[Chapter 3#^1b0173|Eqs. (1.8)]] and [[Chapter 3#^f44962|(1.9)]]

$$
\begin{aligned}
\nabla^{2} F\left(x_{R}\right)=\nabla\left(-R^{\prime}\left(B^{\prime}\right)^{-1} \nabla_{B} f\left(B^{-1}\left(b-R x_{R}\right), x_{R}\right)\right. \\
\left.+\nabla_{R} f\left(B^{-1}\left(b-R x_{R}\right), x_{R}\right)\right)
\end{aligned}
$$

By evaluating this expression at $x_{R}=x_{R}^{*}$, and by partitioning the Hessian $\nabla^{2} f\left(x^{*}\right)$ as

$$
\nabla^{2} f\left(x^{*}\right)=\left(\begin{array}{cc}
\nabla_{B B}^{2} f\left(x^{*}\right) & \nabla_{B R}^{2} f\left(x^{*}\right) \\
\nabla_{R B}^{2} f\left(x^{*}\right) & \nabla_{R R}^{2} f\left(x^{*}\right)
\end{array}\right)
$$

we obtain

$$
\begin{aligned}
\nabla^{2} F\left(x_{R}^{*}\right)= & R^{\prime}\left(B^{\prime}\right)^{-1} \nabla_{B B}^{2} f\left(x^{*}\right) B^{-1} R \\
& -R^{\prime}\left(B^{\prime}\right)^{-1} \nabla_{B R}^{2} f\left(x^{*}\right)-\nabla_{R B}^{2} f\left(x^{*}\right) B^{-1} R+\nabla_{R R}^{2} f\left(x^{*}\right)
\end{aligned}
$$
>[[Note 7 - Second order Hessian derivation for lagrange linear constraints]]

---

This expression together with the positive semidefiniteness of $\nabla^{2} F\left(x_{R}^{*}\right)$ [cf. Eq. (1.14)] and the linearity of the constraints implying that $\nabla^{2} h_{1}\left(x^{*}\right)=0$ , yields for all $d \in \Re^{n-m}$
$$0 \leq d^{\prime} \nabla^{2} F\left(x_{R}^{*}\right) d=y^{\prime} \nabla^{2} f\left(x^{*}\right) y=y^{\prime}\left(\nabla^{2} f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla^{2} h_{i}\left(x^{*}\right)\right) y 
\tag{1.15}$$,
where $y$ is the vector

$$
y=\binom{-B^{-1} R d}{d}
$$

It can be seen that the subspace $V\left(x^{*}\right)$ of feasible variations of Eq. (1.5) is given by

$$
\begin{aligned}
V\left(x^{*}\right) & =\left\{\left(y_{B}, y_{R}\right) \mid B y_{B}+R y_{R}=0\right\} \\
& =\left\{\left(y_{B}, y_{R}\right) \mid y_{B}=-B^{-1} R d, y_{R}=d, d \in \Re^{n-m}\right\}
\end{aligned}
$$

so Eq. (1.15) is equivalent to the second order Lagrange multiplier condition (1.4).

Note that based on the above proof, the Lagrange multiplier conditions of Prop. 3.1.1 are nothing but the "zero gradient" and "positive semidefinite Hessian" conditions for the unconstrained problem (1.8), which is defined over the reduced space of the vector $x_{R}$.

Actually for the case of linear constraints, a stronger version of the Lagrange multiplier theorem can be shown. In particular, the regularity assumption is not needed for the existence of a Lagrange multiplier vector; see Exercise 1.4 and Section 3.4.

We now prove the Lagrange multiplier theorem by generalizing the analysis of the preceding example.

##### Proof of the Lagrange multiplier theorem: 
The proof assumes that $m<n$ : if $m=n$, any vector, including $-\nabla f\left(x^{*}\right)$, can be expressed as a linear combination of the linearly independent vectors $\nabla h_{1}\left(x^{*}\right), \ldots, \nabla h_{m}\left(x^{*}\right)$, thereby proving the theorem. By reordering the coordinates of $x$ if necessary, we partition the vector $x$ as $x=\left(x_{B}, x_{R}\right)$, where the square submatrix $\nabla_{B} h\left(x^{*}\right)$ (the gradient matrix of $h$ with respect to $x_{B}$ ) is invertible. The constraint equation

$$
h\left(x_{B}, x_{R}\right)=0
$$

has the solution $\left(x_{B}^{*}, x_{R}^{*}\right)$, and the implicit function theorem (Prop. A. 25 in Appendix A) can be used to express $x_{B}$ in terms of $x_{R}$ via a unique continuously differentiable function $\phi: S \mapsto \Re^{m}$ defined over a sphere $S$ centered at $x_{R}^{*}$ ( $\phi$ is twice continuously differentiable if $h$ is). In particular, we have $x_{B}^{*}=\phi\left(x_{R}^{*}\right), h\left(\phi\left(x_{R}\right), x_{R}\right)=0$ for all $x_{R} \in S$, and

$$
\nabla \phi\left(x_{R}\right)=-\nabla_{R} h\left(\phi\left(x_{R}\right), x_{R}\right)\left(\nabla_{B} h\left(\phi\left(x_{R}\right), x_{R}\right)\right)^{-1}, \quad \forall x_{R} \in S
$$

---

where $\nabla_{R} h$ is the gradient matrix of $h$ with respect to $x_{R}$.
We now proceed as in the earlier case of linear constraints. We observe that $x_{R}^{*}$ is an unconstrained minimum of the "reduced" cost function

$$
F\left(x_{R}\right)=f\left(\phi\left(x_{R}\right), x_{R}\right)
$$

and we apply the corresponding unconstrained first and second order necessary conditions. The first order Lagrange multiplier condition (1.3) follows by repeating the calculation of Eqs. (1.10)-(1.13) with the definitions

$$
\begin{gathered}
B^{\prime}=\nabla_{B} h\left(x^{*}\right), \quad R^{\prime}=\nabla_{R} h\left(x^{*}\right), \quad A^{\prime}=\nabla h\left(x^{*}\right) \\
\lambda^{*}=-\left(B^{\prime}\right)^{-1} \nabla_{B} f\left(x^{*}\right)
\end{gathered}
$$

The proof of the second order necessary condition (1.4) requires a lengthy calculation, which resembles the one that yielded Eq. (1.15) in the preceding example. In particular, for $d \in \Re^{n-m}$, denote

$$
y=\binom{-B^{-1} R d}{d}
$$

Denote also

$$
H_{i}\left(x_{R}\right)=h_{i}\left(\phi\left(x_{R}\right), x_{R}\right), \quad i=1, \ldots, m
$$

and let $\phi_{i}\left(x_{R}\right)$ be the $i$ th component of $\phi$, so that

$$
\phi\left(x_{R}\right)=\left(\begin{array}{c}
\phi_{1}\left(x_{R}\right) \\
\vdots \\
\phi_{m}\left(x_{R}\right)
\end{array}\right)
$$

Then by differentiating twice the relation

$$
F\left(x_{R}\right)=f\left(\phi\left(x_{R}\right), x_{R}\right)
$$

we obtain by a straightforward calculation [compare with the derivation of Eq. (1.15)]

$$
d^{\prime} \nabla^{2} F\left(x_{R}\right) d=y^{\prime} \nabla^{2} f\left(x^{*}\right) y+d^{\prime}\left(\sum_{j=1}^{m} \nabla^{2} \phi_{j}\left(x_{R}^{*}\right) \frac{\partial f\left(x^{*}\right)}{\partial x_{j}}\right) d
$$

Similarly by twice differentiating the equation $H_{i}\left(x_{R}\right)=0$, we have

$$
0=d^{\prime} \nabla^{2} H_{i}\left(x_{R}\right) d=y^{\prime} \nabla^{2} h_{i}\left(x^{*}\right) y+d^{\prime}\left(\sum_{j=1}^{m} \nabla^{2} \phi_{j}\left(x_{R}^{*}\right) \frac{\partial h_{i}\left(x^{*}\right)}{\partial x_{j}}\right) d
$$

---

By multiplying the above equation with $\lambda_{i}^{*}$, and by adding over all $i$ we obtain

$$
0=\sum_{i=1}^{m} \lambda_{i}^{*} y^{\prime} \nabla^{2} h_{i}\left(x^{*}\right) y+d^{\prime}\left(\sum_{j=1}^{m} \nabla^{2} \mathcal{O}_{j}\left(x_{R}^{*}\right) \sum_{i=1}^{m} \lambda_{i}^{*} \frac{\partial h_{i}\left(x^{*}\right)}{\partial x_{j}}\right) d
$$

Adding Eqs. (1.16) and (1.17), and using the relations $d^{\prime} \nabla^{2} F\left(x_{R}\right) d \geq 0$ and

$$
\frac{\partial f\left(x^{*}\right)}{\partial x_{j}}+\sum_{i=1}^{m} \lambda_{i}^{*} \frac{\partial h_{i}\left(x^{*}\right)}{\partial x_{j}}=0, \quad j=1, \ldots, m
$$

we obtain

$$
0 \leq y^{\prime}\left(\nabla^{2} f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla^{2} h_{i}\left(x^{*}\right)\right) y
$$

for all $y$ of the form defined above. As shown in Example 1.1, $y$ belongs to the subspace $U\left(x^{*}\right)$ if and only if $y$ has this form, thus proving the second order Lagrange multiplier condition (1.4). Q.E.D.

# 3.1.2 The Penalty Approach 

The line of analysis here is based on approximating the original constrained problem by an unconstrained optimization problem that involves a penalty for violation of the constraints, and then applying the unconstrained optimality conditions to the latter problem. In particular, for $k=1,2, \ldots$, we introduce the cost function

$$
F^{k}(x)=f(x)+\frac{k}{2}\|h(x)\|^{2}+\frac{\alpha}{2}\left\|x-x^{*}\right\|^{2}
$$

where $x^{*}$ is the local minimum of the constrained problem and $\alpha$ is some positive scalar. The term $(k / 2)\|h(x)\|^{2}$ imposes a penalty for violating the constraint $h(x)=0$, while the term $(\alpha / 2)\left\|x-x^{*}\right\|^{2}$ is introduced for technical proof-related reasons |to ensure that $x^{*}$ is a strict local minimum of the function $f(x)+(\alpha / 2)\left\|x-x^{*}\right\|^{2}$ subject to $\left.h(x)=0\right\|$.

Since $x^{*}$ is a local minimum, we can select $\epsilon>0$ such that $f\left(x^{*}\right) \leq$ $f(x)$ for all feasible $x$ in the closed sphere

$$
S=\left\{x \mid\left\|x-x^{*}\right\| \leq \epsilon\right\}
$$

Let $x^{k}$ be an optimal solution of the problem

$$
\begin{aligned}
& \text { minimize } F^{k}(x) \\
& \text { subject to } x \in S
\end{aligned}
$$

---

[An optimal solution exists because of Weierstrass' theorem (Prop. A. 8 in Appendix A), since $S$ is compact.] We will show that the sequence $\left\{x^{k}\right\}$ converges to $x^{*}$.

We have

$$
F^{k}\left(x^{k}\right)=f\left(x^{k}\right)+\frac{k}{2}\left\|h\left(x^{k}\right)\right\|^{2}+\frac{\alpha}{2}\left\|x^{k}-x^{*}\right\|^{2} \leq F^{k}\left(x^{*}\right)=f\left(x^{*}\right)
$$

and since $f\left(x^{k}\right)$ and $\left\|x^{k}-x^{*}\right\|^{2}$ are bounded over $S$, we obtain

$$
\lim _{k \rightarrow \infty}\left\|h\left(x^{k}\right)\right\|=0
$$

Therefore, for every limit point $\bar{x}$ of $\left\{x^{k}\right\}$ we have $h(\bar{x})=0$. Furthermore, Eq. (1.19) yields $f\left(x^{k}\right)+(\alpha / 2)\left\|x^{k}-x^{*}\right\|^{2} \leq f\left(x^{*}\right)$ for all $k$, so by taking the limit as $k \rightarrow \infty$, we obtain

$$
f(\bar{x})+\frac{\alpha}{2}\left\|\bar{x}-x^{*}\right\|^{2} \leq f\left(x^{*}\right)
$$

Since $\bar{x} \in S$ and $\bar{x}$ is feasible, we have $f\left(x^{*}\right) \leq f(\bar{x})$, which when combined with the preceding inequality yields $\left\|\bar{x}-x^{*}\right\|=0$ so that $\bar{x}=x^{*}$. Thus the sequence $\left\{x^{k}\right\}$ converges to $x^{*}$, and it follows that $x^{k}$ is an interior point of the closed sphere $S$ for sufficiently large $k$. Therefore $x^{k}$ is an unconstrained local minimum of $F^{k}(x)$ for sufficiently large $k$. We will now prove the Lagrange multiplier theorem by working with the corresponding unconstrained necessary optimality conditions.

Proof of the Lagrange multiplier theorem: From the first order necessary condition, we have

$$
0=\nabla F^{k}\left(x^{k}\right)=\nabla f\left(x^{k}\right)+k \nabla h\left(x^{k}\right) h\left(x^{k}\right)+\alpha\left(x^{k}-x^{*}\right)
$$

Since $\nabla h\left(x^{*}\right)$ has rank $m$, the same is true for $\nabla h\left(x^{k}\right)$ if $k$ is sufficiently large. For such $k, \nabla h\left(x^{k}\right)^{\prime} \nabla h\left(x^{k}\right)$ is invertible (Prop. A. 20 in Appendix A), and by premultiplying Eq. (1.20) with $\left(\nabla h\left(x^{k}\right)^{\prime} \nabla h\left(x^{k}\right)\right)^{-1} \nabla h\left(x^{k}\right)^{\prime}$, we obtain

$$
k h\left(x^{k}\right)=-\left(\nabla h\left(x^{k}\right)^{\prime} \nabla h\left(x^{k}\right)\right)^{-1} \nabla h\left(x^{k}\right)^{\prime}\left(\nabla f\left(x^{k}\right)+\alpha\left(x^{k}-x^{*}\right)\right)
$$

By taking the limit as $k \rightarrow \infty$ and $x^{k} \rightarrow x^{*}$, we see that $\left\{k h\left(x^{k}\right)\right\}$ converges to the vector

$$
\lambda^{*}=-\left(\nabla h\left(x^{*}\right)^{\prime} \nabla h\left(x^{*}\right)\right)^{-1} \nabla h\left(x^{*}\right)^{\prime} \nabla f\left(x^{*}\right)
$$

By taking the limit as $k \rightarrow \infty$ in Eq. (1.20), we obtain

$$
\nabla f\left(x^{*}\right)+\nabla h\left(x^{*}\right) \lambda^{*}=0
$$

---

proving the first order Lagrange multiplier condition (1.3).
By using the second order unconstrained optimality condition for problem (1.18), we see that the matrix

$$
\nabla^{2} F^{k}\left(x^{k}\right)=\nabla^{2} f\left(x^{k}\right)+k \nabla h\left(x^{k}\right) \nabla h\left(x^{k}\right)^{\prime}+k \sum_{i=1}^{m} h_{i}\left(x^{k}\right) \nabla^{2} h_{i}\left(x^{k}\right)+\alpha I
$$

is positive semidefinite, for all $k$ and $\alpha>0$. Fix any $y \in V\left(x^{*}\right)$ [that is, $\nabla h\left(x^{*}\right)^{\prime} y=0$ ], and let $y^{k}$ be the projection of $y$ on the nullspace of $\nabla h\left(x^{k}\right)^{\prime}$, that is,

$$
y^{k}=y-\nabla h\left(x^{k}\right)\left(\nabla h\left(x^{k}\right)^{\prime} \nabla h\left(x^{k}\right)\right)^{-1} \nabla h\left(x^{k}\right)^{\prime} y
$$

[cf. Eq. (1.24) in Section 2.1]. Since $\nabla h\left(x^{k}\right)^{\prime} y^{k}=0$ and $\nabla^{2} F^{k}\left(x^{k}\right)$ is positive semidefinite, we have
$0 \leq y^{k^{\prime}} \nabla^{2} F^{k}\left(x^{k}\right) y^{k}=y^{k^{\prime}}\left(\nabla^{2} f\left(x^{k}\right)+k \sum_{i=1}^{m} h_{i}\left(x^{k}\right) \nabla^{2} h_{i}\left(x^{k}\right)\right) y^{k}+\alpha\left\|y^{k}\right\|^{2}$.
Since $k h_{i}\left(x^{k}\right) \rightarrow \lambda_{i}^{*}$, and from Eq. (1.21), together with the facts $x^{k} \rightarrow x^{*}$ and $\nabla h\left(x^{*}\right)^{\prime} y=0$, we have $y^{k} \rightarrow y$, we obtain

$$
0 \leq y^{\prime}\left(\nabla^{2} f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla^{2} h_{i}\left(x^{*}\right)\right) y+\alpha\|y\|^{2}, \quad \forall y \in V\left(x^{*}\right)
$$

Since $\alpha$ can be taken arbitrarily close to 0 , we obtain

$$
0 \leq y^{\prime}\left(\nabla^{2} f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla^{2} h_{i}\left(x^{*}\right)\right) y, \quad \forall y \in V\left(x^{*}\right)
$$

which is the second order Lagrange multiplier condition. The proof of the Lagrange multiplier theorem is thus complete. Q.E.D.

We now illustrate what can go wrong when the regularity assumption in the Lagrange multiplier theorem is violated.

# Example 1.2 (A Problem with no Lagrange Multipliers) 

Consider the problem

$$
\begin{aligned}
\operatorname{minimize} & f(x)=x_{1}+x_{2} \\
\text { subject to } & h_{1}(x)=\left(x_{1}-1\right)^{2}+x_{2}^{2}-1=0 \\
& h_{2}(x)=\left(x_{1}-2\right)^{2}+x_{2}^{2}-4=0
\end{aligned}
$$

---

The geometry of this problem is illustrated in Fig. 3.1.2. It can be seen that at the local minimum $x^{*}=(0,0)$, the cost gradient $\nabla f\left(x^{*}\right)=(1,1)$ cannot be expressed as a linear combination of the constraint gradients $\nabla h_{1}\left(x^{*}\right)=$ $(-2,0)$ and $\nabla h_{2}\left(x^{*}\right)=(-4,0)$. Thus the Lagrange multiplier condition

$$
\nabla f\left(x^{*}\right)+\lambda_{1}^{*} \nabla h_{1}\left(x^{*}\right)+\lambda_{2}^{*} \nabla h_{2}\left(x^{*}\right)=0
$$

cannot hold for any $\lambda_{1}^{*}$ and $\lambda_{2}^{*}$.
The difficulty here is that the subspace of first order feasible variations [cf. Eq. (1.5)], which is $\left\{y \mid y_{1}=0, y_{2} \in \Re\right\}$, has larger dimension than the true set of feasible variations $\{y \mid y=0\}$. The optimality of $x^{*}$ implies that $\nabla f\left(x^{*}\right)$ is orthogonal to the true set of feasible variations, but for a Lagrange multiplier to exist, $\nabla f\left(x^{*}\right)$ must be orthogonal to the subspace of first order feasible variations. If the constraint gradients $\nabla h_{1}\left(x^{*}\right)$ and $\nabla h_{2}\left(x^{*}\right)$ were linearly independent, the set of feasible variations and the subspace of first order feasible variations would have the same dimension and this difficulty would not have occurred.
![[ch3_and_4_p12_img2.jpeg]]

Figure 3.1.2. Illustration of how Lagrange multipliers may not exist. The problem here is

$$
\begin{aligned}
\operatorname{minimize} & f(x)=x_{1}+x_{2} \\
\text { subject to } & h_{1}(x)=\left(x_{1}-1\right)^{2}+x_{2}^{2}-1=0 \\
& h_{2}(x)=\left(x_{1}-2\right)^{2}+x_{2}^{2}-4=0
\end{aligned}
$$

with a local minimum $x^{*}=(0,0)$. The cost gradient cannot be expressed as a linear combination of the constraint gradients, so there are no Lagrange multipliers.

---

# 3.1.3 The Lagrangian Function 

Sometimes it is convenient to write our necessary conditions in terms of the Lagrangian function $L: \Re^{n+m} \mapsto \Re$ defined by

$$
L(x, \lambda)=f(x)+\sum_{i=1}^{m} \lambda_{i} h_{i}(x)
$$

Then, if $x^{*}$ is a local minimum which is regular, the Lagrange multiplier conditions of Prop. 3.1.1 together with the equation $h\left(x^{*}\right)=0$ are written compactly as

$$
\begin{gathered}
\nabla_{x} L\left(x^{*}, \lambda^{*}\right)=0, \quad \nabla_{\lambda} L\left(x^{*}, \lambda^{*}\right)=0 \\
y^{\prime} \nabla_{x x}^{2} L\left(x^{*}, \lambda^{*}\right) y \geq 0, \quad \text { for all } y \in V\left(x^{*}\right)
\end{gathered}
$$

The first order necessary conditions (1.23) represent a system of $n+m$ equations with $n+m$ unknowns the coordinates of $x^{*}$ and $\lambda^{*}$. Every local minimum $x^{*}$ which is regular, together with its associated Lagrange multiplier vector, will be a solution of this system. However, a solution of the system need not correspond to a local minimum, as our experience with unconstrained problems indicates: it could be for example a local maximum.

## Example 1.3

Consider the problem

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2}\left(x_{1}^{2}+x_{2}^{2}+x_{3}^{2}\right) \\
& \text { subject to } x_{1}+x_{2}+x_{3}=3
\end{aligned}
$$

Then the first order necessary conditions (1.23) yield

$$
\begin{aligned}
x_{1}^{*}+\lambda^{*} & =0 \\
x_{2}^{*}+\lambda^{*} & =0 \\
x_{3}^{*}+\lambda^{*} & =0 \\
x_{1}^{*}+x_{2}^{*}+x_{3}^{*} & =3
\end{aligned}
$$

This is a system of four equations and four unknowns. It has the unique solution

$$
x_{1}^{*}=x_{2}^{*}=x_{3}^{*}=1, \quad \lambda^{*}=-1
$$

The constraint gradient here is $(1,1,1)$, so all feasible vectors are regular. Therefore, $x^{*}=(1,1,1)$ is the unique candidate for a local minimum. Furthermore, since $\nabla_{x x}^{2} L\left(x^{*}, \lambda^{*}\right)$ is the identity matrix for this problem, the second order necessary condition (1.24) is satisfied.

---

If instead we consider the problem

$$
\begin{array}{ll}
\operatorname{minimize} & -\frac{1}{2}\left(x_{1}^{2}+x_{2}^{2}+x_{3}^{2}\right) \\
\text { subject to } & x_{1}+x_{2}+x_{3}=3
\end{array}
$$

then the first order condition (1.23) yields $x^{*}=(1,1,1)$ and $\lambda^{*}=1$. However, the second order condition (1.24) is not satisfied and since every feasible vector is also regular, we conclude that there is no local minimum.

In fact, by using the convexity (or concavity) of the cost function of problem (1.25) [or problem (1.26), respectively] and the convexity of the constraint set, it can be seen using Prop. 2.1.1 that $x^{*}$ is a global minimum (maximum, respectively).

# Example 1.4 (Diffraction Law in Optics) 

Consider a smooth curve on the 2-dimensional plane described by the equation

$$
h(x)=0
$$

where $h: \Re^{2} \mapsto \Re$ is continuously differentiable. Imagine that the curve separates the plane into two regions and that the velocity of light is different in each region. Let $y$ and $z$ be two points that lie on opposite sides of the curve as shown in part (a) of Fig. 3.1.3. Suppose that light going from $y$ to $z$ follows a ray from $y$ to some point $x^{*}$ of the curve with velocity $v_{y}$, and then follows a ray from $x^{*}$ to $z$ with velocity $v_{z}$. Assuming that light follows a path of minimum travel time, we will show that $x^{*}$ is characterized by the following diffraction law:

$$
\frac{\sin \phi_{y}}{v_{y}}=\frac{\sin \phi_{z}}{v_{z}}
$$

where $\phi_{y}$ and $\phi_{z}$ are the angles shown in Fig. 3.1.3.
Indeed, the travel time can be expressed as

$$
T\left(x^{*}\right)=\frac{\left\|y-x^{*}\right\|}{v_{y}}+\frac{\left\|z-x^{*}\right\|}{v_{z}}
$$

and the Lagrange multiplier theorem states that the gradient

$$
\nabla T\left(x^{*}\right)=\frac{x^{*}-y}{v_{y}\left\|x^{*}-y\right\|}+\frac{x^{*}-z}{v_{z}\left\|x^{*}-z\right\|}
$$

is a scalar multiple of $\nabla h\left(x^{*}\right)$. [For this we need that $x^{*}$ is regular, that is, $\nabla h\left(x^{*}\right) \neq 0$, which we assume.] To interpret this condition, consider the vectors

$$
\bar{y}=x^{*}+\frac{y-x^{*}}{v_{y}\left\|y-x^{*}\right\|}, \quad \bar{z}=x^{*}+\frac{z-x^{*}}{v_{z}\left\|z-x^{*}\right\|}
$$

and the parallelogram with sides $\bar{y}-x^{*}$ and $\bar{z}-x^{*}$, which by Eq. (1.28), has $-\nabla T\left(x^{*}\right)$ as one of its diagonals. From part (b) of Fig. 3.1.3, we see that the Lagrange multiplier condition states that the diagonal $-\nabla T\left(x^{*}\right)$ of the parallelogram is normal to the curve at $x^{*}$. From Euclidean geometry it follows that the distances of $\bar{y}$ and $\bar{z}$ from the vertical diagonal are equal. Since $\left\|\bar{y}-x^{*}\right\|=1 / v_{y}$ and $\left\|\bar{z}-x^{*}\right\|=1 / v_{z}$, this is equivalent to the diffraction law (1.27).

---

![[ch3_and_4_p15_img3.jpeg]]

Figure 3.1.3. Diffraction law; cf. Example 1.4. The angles $\phi_{y}$ and $\phi_{z}$ are formed by the normal line to the curve at $x^{*}$ [which is $\nabla h\left(x^{*}\right)$ ], and the vectors $y-x^{*}$ and $z-x^{*}$, respectively.

# Example 1.5 (Optimal Portfolio Selection) 

Consider an investor who wants to allocate one unit of wealth among $n$ assets offering random rates of return $e_{1}, \ldots, e_{n}$, respectively. The means $\bar{e}_{i}=$ $E\left\{e_{i}\right\}, i=1, \ldots, n$, and the covariance matrix

$$
Q=\left(\begin{array}{ccc}
E\left\{\left(e_{1}-\bar{e}_{1}\right)^{2}\right\} & \ldots & E\left\{\left(e_{1}-\bar{e}_{1}\right)\left(e_{n}-\bar{e}_{n}\right)\right\} \\
\vdots & \ldots & \vdots \\
E\left\{\left(e_{n}-\bar{e}_{n}\right)\left(e_{1}-\bar{e}_{1}\right)\right\} & \ldots & E\left\{\left(e_{n}-\bar{e}_{n}\right)^{2}\right\}
\end{array}\right)
$$

are known, and we assume that $Q$ is invertible. If $x_{i}$ is the amount invested in the ith asset, the mean and the variance of the return of the investment $y=\sum_{i=1}^{n} e_{i} x_{i}$ are

$$
\bar{y}=E\{y\}=\sum_{i=1}^{n} \bar{e}_{i} x_{i}
$$

and

$$
\begin{aligned}
\sigma^{2} & =E\left\{(y-\bar{y})^{2}\right\}=E\left\{\left(\sum_{i=1}^{n}\left(e_{i}-\bar{e}_{i}\right) x_{i}\right)^{2}\right\} \\
& =\sum_{i=1}^{n} \sum_{j=1}^{n} E\left\{\left(e_{i}-\bar{e}_{i}\right)\left(e_{j}-\bar{e}_{j}\right)\right\} x_{i} x_{j}=x^{\prime} Q x
\end{aligned}
$$

The investor's problem is to find the portfolio $x=\left(x_{1}, \ldots x_{n}\right)$ that minimizes the variance $E\left\{(y-\bar{y})^{2}\right\}$ to achieve a given level of mean return $E\{y\}$, say $E\{y\}=m$. Thus the problem is

$$
\begin{aligned}
& \text { minimize } x^{\prime} Q x \\
& \text { subject to } \sum_{i=1}^{n} x_{i}=1, \quad \sum_{i=1}^{n} \bar{e}_{i} x_{i}=m
\end{aligned}
$$

---

We want to see how the solution varies with $m$. Note that the solution of the problem "scales" with the amount of wealth available for investment. In particular, it can be seen that if $w$ is the amount available for investment and $w m$ is the desired mean, then the optimal portfolio is $w x^{*}$, where $x^{*}$ is the optimal portfolio corresponding to one unit of wealth and a desired mean of $m$.

Let $\lambda_{1}$ and $\lambda_{2}$ be Lagrange multipliers for the constraints $\sum_{i=1}^{n} x_{i}=1$ and $\sum_{i=1}^{n} \bar{e}_{i} x_{i}=m$, respectively. The first order optimality condition is

$$
2 Q x^{*}+\lambda_{1} u+\lambda_{2} \bar{e}=0
$$

where $u=(1, \ldots, 1)^{\prime}$ and $\bar{e}=\left(\bar{e}_{1}, \ldots, \bar{e}_{n}\right)^{\prime}$. Equivalently,

$$
x^{*}=-\frac{1}{2} Q^{-1} u \lambda_{1}-\frac{1}{2} Q^{-1} \bar{e} \lambda_{2}
$$

and by substitution in the constraints $u^{\prime} x^{*}=1$ and $\bar{e}^{\prime} x^{*}=m$, we obtain

$$
\begin{aligned}
& 1=u^{\prime} x^{*}=-\frac{1}{2} u^{\prime} Q^{-1} u \lambda_{1}-\frac{1}{2} u^{\prime} Q^{-1} \bar{e} \lambda_{2} \\
& m=\bar{e}^{\prime} x^{*}=-\frac{1}{2} \bar{e}^{\prime} Q^{-1} u \lambda_{1}-\frac{1}{2} \bar{e}^{\prime} Q^{-1} \bar{e} \lambda_{2}
\end{aligned}
$$

By solving these equations for $\lambda_{1}$ and $\lambda_{2}$ we obtain

$$
\begin{aligned}
& \lambda_{1}=\xi_{1}+\zeta_{1} m \\
& \lambda_{2}=\xi_{m}+\zeta_{m} m
\end{aligned}
$$

for some scalars $\xi_{1}, \zeta_{1}, \xi_{m}, \zeta_{m}$. Thus substitution in Eq. (1.29) yields

$$
x^{*}=m v+w
$$

for some vectors $v$ and $w$ that depend on $Q$ and $\bar{e}$. The corresponding variance of return is

$$
\sigma^{2}=(m v+w)^{\prime} Q(m v+w)=(\alpha m+\beta)^{2}+\gamma
$$

where $\alpha, \beta$, and $\gamma$ are some scalars that depend on $Q$ and $\bar{e}$.
Equation (1.31) specifies the locus of pairs $(\sigma, m)$ that are achievable by optimal portfolio selection. Suppose now that one of the assets is riskless, that is, its return is fixed at some known value $\bar{e}_{f}$ and the variance of its return is zero. Then $\gamma$ must be zero since when $m=\bar{e}_{f}$, the variance is minimized (set to zero) by investing exclusively in the riskless asset. In this case, from Eq. (1.31) we obtain

$$
\sigma=|\alpha m+\beta|
$$

Thus, assuming $\alpha>0$, portfolios of interest correspond to mean-variance pairs of the form

$$
\sigma=\alpha m+\beta, \quad m \geq \bar{e}_{f}
$$

---

as shown in Fig. 3.1.4. This is known as the efficient frontier. To characterize the efficient frontier it is sufficient to determine one more point on it. Every other point is determined as a linear combination of this point and the point $(0, \bar{\epsilon})$. Furthermore, by Eq. (1.30), every pair $(\sigma, m)$ on the efficient frontier can be achieved by an appropriate linear combination of two portfolios, one consisting of exclusive investment on the riskless asset, and the other corresponding to some pair $(\sigma, m)$ on the efficient frontier with $\sigma \neq 0$. (Economic arguments can be used to construct a second portfolio on the efficient frontier, leading to the celebrated CAPM model of finance theory; see e.g. [Hul,88]. Thus, while different investors with different preferences towards risk may prefer different points on the efficient frontier, their preferences can be obtained by appropriate mixtures of just two portfolios.
![[ch3_and_4_p17_img4.jpeg]]

Figure 3.1.4. Efficient frontier for optimal portfolio selection. Every portfolio corresponding to a mean-variance pair on the efficient frontier can be achieved as a linear combination of two portfolios corresponding to two pairs on the efficient frontier.

# E X E R C ISES 

## 1.1

Use the Lagrange multiplier theorem to solve the following problems:
(a) $f(x)=\|x\|^{2}, h(x)=\sum_{i=1}^{n} x_{i}-1$.
(b) $f(x)=\sum_{i=1}^{n} x_{i}, h(x)=\|x\|^{2}-1$.
(c) $f(x)=\|x\|^{2}, h(x)=x^{\prime} Q x-1$, where $Q$ is positive definite.

## 1.2

Among all pyramids whose base is a given triangle and whose height is given, find the one that has maximal surface area (sum of areas of the triangular sides). Answer: The optimal vertex projects on the point in the triangle that is at equal distance from all the sides.

---

# 1.3 (Fermat's Principle in Optics) 

Consider a smooth curve on the plane described by the equation

$$
h(x)=0
$$

where $h: \Re^{2} \mapsto \Re$ is continuously differentiable. Let $y$ and $z$ be two points that lie in relation to the curve as shown in Fig. 3.1.5. Show that if a point $x^{*}$ minimizes the sum of the Euclidean distances

$$
\|y-x\|+\|z-x\|
$$

over all points $x$ on the curve, then the angles $\phi_{y}, \phi_{z}$ shown in Fig. 3.1.6 must be equal.
![[ch3_and_4_p18_img5.jpeg]]

Figure 3.1.5. Fermat's principle; cf. Exercise 1.3. The angles $\phi_{y}, \phi_{z}$ must be equal.

## 1.4

Show that if the constraints are linear, the regularity assumption is not needed in order for the Lagrange multiplier conditions to hold, except that the Lagrange multiplier vector need not be unique. Hint: Discard the redundant equality constraints, and assign to them zero Lagrange multipliers.

## 1.5

Consider a symmetric $n \times n$ matrix $Q$. Define

$$
\lambda_{1}=\min _{\|x\|^{2}=1} x^{\prime} Q x, \quad e_{1}=\arg \min _{\|x\|^{2}=1} x^{\prime} Q x
$$

and for $k=0, \ldots, n-1$,

$$
\lambda_{k+1}=\min _{\substack{\|x\|^{2}=1 \\ e_{1}^{\prime} x=0, i=1, \ldots, k}} x^{\prime} Q x, \quad e_{k+1}=\arg \min _{\substack{\|x\|^{2}=1 \\ e_{1}^{\prime} x=0, i=1, \ldots, k}} x^{\prime} Q x
$$

---

(a) Show that

$$
\lambda_{1} \leq \lambda_{2} \leq \cdots \leq \lambda_{n}
$$

(b) Show that the vectors $e_{1}, \ldots, e_{n}$ are linearly independent.
(c) Interpret $\lambda_{1}, \ldots, \lambda_{n}$ as Lagrange multipliers, and show that $\lambda_{1}, \ldots, \lambda_{n}$ are eigenvalues of $Q$, while $e_{1}, \ldots, e_{n}$ are corresponding eigenvectors.

# 1.6 

Solve the problem

$$
\begin{aligned}
& \operatorname{minimize} \sum_{i=1}^{n} x_{i} \\
& \text { subject to } x_{1} x_{2} \cdots x_{n}=1, \quad x_{i}>0, \quad i=1, \ldots, n
\end{aligned}
$$

Establish the arithmetic-geometric mean inequality

$$
\left(x_{1} x_{2} \cdots x_{n}\right)^{1 / n} \leq \frac{\sum_{i=1}^{n} x_{i}}{n}
$$

for a set of positive numbers $x_{i}, i=1, \ldots, n$.

## 1.7

Show that the minimum of the three-dimensional problem

$$
\begin{aligned}
& \text { minimize } x+y+z \\
& \text { subject to } x y z=1
\end{aligned}
$$

is attained at four distinct vectors.

### 1.2 SUFFICIENT CONDITIONS AND SENSITIVITY ANALYSIS

As shown by Example 1.3 of the preceding section, the first order necessary condition may be satisfied by both local minima and local maxima (and possibly other vectors). The second order necessary condition is useful in narrowing down the field of candidates for local minima. To guarantee that a given vector is a local minimum, we need sufficient conditions for optimality.

One way to obtain sufficient optimality conditions is to use the elimination approach of the previous section, and convert the problem to the

---

unconstrained problem $\min _{x_{R} \in \Re^{n-m}} F\left(x_{R}\right)$, where $F$ is the "reduced" cost function. One may then write the sufficient conditions $\nabla F\left(x_{R}^{*}\right)=0$ and $\nabla^{2} F\left(x_{R}^{*}\right)$ : positive definite, and translate them to the constrained context. The drawback is that regularity of $x^{*}$ must be assumed in order for the elimination approach to be valid. It turns out that the same conditions can be obtained directly without the regularity assumption, as shown in the following proposition.

Proposition 3.2.1: (Second Order Sufficiency Conditions) Assume that $f$ and $h$ are twice continuously differentiable, and let $x^{*} \in$ $\Re^{n}$ and $\lambda^{*} \in \Re^{m}$ satisfy

$$
\begin{gathered}
\nabla_{x} L\left(x^{*}, \lambda^{*}\right)=0, \quad \nabla_{\lambda} L\left(x^{*}, \lambda^{*}\right)=0 \\
y^{\prime} \nabla_{x x}^{2} L\left(x^{*}, \lambda^{*}\right) y>0, \quad \text { for all } y \neq 0 \text { with } \nabla h\left(x^{*}\right)^{\prime} y=0
\end{gathered}
$$

Then $x^{*}$ is a strict local minimum of $f$ subject to $h(x)=0$. In fact, there exist scalars $\gamma>0$ and $\epsilon>0$ such that
$f(x) \geq f\left(x^{*}\right)+\frac{\gamma}{2}\left\|x-x^{*}\right\|^{2}, \quad \forall x$ with $h(x)=0$ and $\left\|x-x^{*}\right\|<\epsilon$.

Proof: We first note that $x^{*}$ is feasible, since, by Eq. (2.1), $\nabla_{\lambda} L\left(x^{*}, \lambda^{*}\right)=$ $h\left(x^{*}\right)=0$. Assume that the conclusion of the proposition does not hold, so that there is a sequence $\left\{x^{k}\right\}$ such that $x^{k} \rightarrow x^{*}$, and for all $k, x^{k} \neq x^{*}$, $h\left(x^{k}\right)=0$, and $f\left(x^{k}\right)<f\left(x^{*}\right)+(1 / k)\left\|x^{k}-x^{*}\right\|^{2}$. Let us write $x^{k}=$ $x^{*}+\delta^{k} y^{k}$, where

$$
\delta^{k}=\left\|x^{k}-x^{*}\right\|, \quad y^{k}=\frac{x^{k}-x^{*}}{\left\|x^{k}-x^{*}\right\|}
$$

The sequence $\left\{y^{k}\right\}$ is bounded, so it must have a subsequence converging to some $y$ with $\|y\|=1$. Without loss of generality, we assume that the whole sequence $\left\{y^{k}\right\}$ converges to $y$. By taking the limit as $\delta^{k} \rightarrow 0$ in the relation

$$
0=\frac{h_{i}\left(x^{k}\right)-h_{i}\left(x^{*}\right)}{\delta^{k}}=\frac{h_{i}\left(x^{*}+\delta^{k} y^{k}\right)-h_{i}\left(x^{*}\right)}{\delta^{k}}=\nabla h_{i}\left(x^{*}\right)^{\prime} y+\frac{o\left(\delta^{k}\right)}{\delta^{k}}
$$

we see that $\nabla h\left(x^{*}\right)^{\prime} y=0$.
We will now show that $y^{\prime} \nabla_{x x}^{2} L\left(x^{*}, \lambda^{*}\right) y \leq 0$, thus coming to a contradiction [cf. Eq. (2.2)]. Since $x_{k}=x^{*}+\delta^{k} y^{k}$, by the mean value theorem [Prop. A.23(b) in Appendix A], we have

$$
0=h_{i}\left(x^{k}\right)-h_{i}\left(x^{*}\right)=\delta^{k} \nabla h_{i}\left(x^{*}\right)^{\prime} y^{k}+\frac{\left(\delta^{k}\right)^{2}}{2} y^{k^{\prime}} \nabla^{2} h_{i}\left(\tilde{\xi}_{i}^{k}\right) y^{k}
$$

---

$$
\frac{1}{k}\left\|x^{k}-x^{*}\right\|^{2}>f\left(x^{k}\right)-f\left(x^{*}\right)=\delta^{k} \nabla f\left(x^{*}\right)^{\prime} y^{k}+\frac{\left(\delta^{k}\right)^{2}}{2} y^{k^{\prime}} \nabla^{2} f\left(\tilde{\xi}^{k}\right) y^{k}
$$

where each of the vectors $\tilde{\xi}_{i}^{k}$ and $\tilde{\xi}^{k}$ lies on the line segment joining $x^{*}$ and $x^{k}$. Multiplying Eq. (2.3) by $\lambda_{i}^{*}$ and adding Eq. (2.4) to it, we obtain

$$
\begin{aligned}
\frac{1}{k}\left\|x^{k}-x^{*}\right\|^{2}> & \delta^{k}\left(\nabla f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla h_{i}\left(x^{*}\right)\right)^{\prime} y^{k} \\
& +\frac{\left(\delta^{k}\right)^{2}}{2} y^{k^{\prime}}\left(\nabla^{2} f\left(\tilde{\xi}^{k}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla^{2} h_{i}\left(\tilde{\xi}_{i}^{k}\right)\right) y^{k}
\end{aligned}
$$

Since $\delta^{k}=\left\|x^{k}-x^{*}\right\|$ and $\nabla f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla h_{i}\left(x^{*}\right)=0$, we obtain

$$
\frac{2}{k}>y^{k^{\prime}}\left(\nabla^{2} f\left(\tilde{\xi}^{k}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla^{2} h_{i}\left(\tilde{\xi}_{i}^{k}\right)\right) y^{k}
$$

By taking the limit as $k \rightarrow \infty$,

$$
0 \geq y^{\prime}\left(\nabla^{2} f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla^{2} h_{i}\left(x^{*}\right)\right) y
$$

which contradicts the hypothesis (2.2). Q.E.D.

# Example 2.1 

Consider the problem

$$
\begin{aligned}
& \text { minimize }-\left(x_{1} x_{2}+x_{2} x_{3}+x_{1} x_{3}\right) \\
& \text { subject to } x_{1}+x_{2}+x_{3}=3
\end{aligned}
$$

Note that if $x_{1}, x_{2}$, and $x_{3}$ represent the length, width, and height of a parallelepiped $P$, respectively, the problem can be interpreted as maximizing the surface area of $P$, subject to the sum of the edge lengths of $P$ being fixed. The first order necessary conditions are

$$
\begin{aligned}
-x_{2}^{*}-x_{3}^{*}+\lambda^{*}=0 \\
-x_{1}^{*}-x_{3}^{*}+\lambda^{*}=0 \\
-x_{1}^{*}-x_{2}^{*}+\lambda^{*}=0 \\
x_{1}^{*}+x_{2}^{*}+x_{3}^{*}=3
\end{aligned}
$$

---

which have the unique solution $x_{1}^{*}=x_{2}^{*}=x_{3}^{*}=1, \lambda^{*}=2$. The Hessian of the Lagrangian is

$$
\nabla_{x x}^{2} L\left(x^{*}, \lambda^{*}\right)=\left(\begin{array}{ccc}
0 & -1 & -1 \\
-1 & 0 & -1 \\
-1 & -1 & 0
\end{array}\right)
$$

We have for all $y \in V=\left\{y \mid \nabla h\left(x^{*}\right)^{\prime} y=0\right\}=\left\{y \mid y_{1}+y_{2}+y_{3}=0\right\}$ with $y \neq 0$,

$$
\begin{aligned}
y^{\prime} \nabla_{x x}^{2} L\left(x^{*}, \lambda^{*}\right) y & =-y_{1}\left(y_{2}+y_{3}\right)-y_{2}\left(y_{1}+y_{3}\right)-y_{3}\left(y_{1}+y_{2}\right) \\
& =y_{1}^{2}+y_{2}^{2}+y_{3}^{2}>0
\end{aligned}
$$

Hence, the sufficient conditions of Prop. 3.2.1 are satisfied and $x^{*}$ is a strict local minimum.

# 3.2.1 Proofs Using the Augmented Lagrangian 

It is instructive to provide an alternative proof of the sufficiency conditions, using concepts that will be used later in Section 4.2 as the basis for some important algorithms.

We first show a useful lemma. Consider a symmetric matrix $P$ and a positive semidefinite matrix $Q$. For each $x$ not in the nullspace of $Q$, we have $x^{\prime} Q x>0$, so $x^{\prime}(P+c Q) x>0$ for sufficiently large scalars $c$. Thus, if we assume that $x^{\prime} P x>0$ for all $x \neq 0$ in the nullspace of $Q$, then for every $x \neq 0$, we have $x^{\prime}(P+c Q) x>0$, provided that $c$ is sufficiently large. The essence of the lemma is that we can choose a threshold value of $c$ that works for all $x$, so that $P+c Q$ is positive definite for all $c$ greater than that threshold.

Lemma 3.2.1: Let $P$ and $Q$ be two symmetric matrices. Assume that $Q$ is positive semidefinite and $P$ is positive definite on the nullspace of $Q$, that is, $x^{\prime} P x>0$ for all $x \neq 0$ with $x^{\prime} Q x=0$. Then there exists a scalar $\bar{c}$ such that

$$
P+c Q: \text { positive definite, } \quad \forall c>\bar{c}
$$

Proof: Assume the contrary. Then for every integer $k$, there exists a vector $x^{k}$ with $\left\|x^{k}\right\|=1$ such that

$$
x^{k^{\prime}} P x^{k}+k x^{k^{\prime}} Q x^{k} \leq 0
$$

---

Since $\left\{x^{k}\right\}$ is bounded, there is a subsequence $\left\{x^{k}\right\}_{k \in K}$ converging to some $\bar{x}$ [Prop. A.5(c) in Appendix A], and since $\left\|x^{k}\right\|=1$ for all $k$, we have $\|\bar{x}\|=1$. Taking the limit superior in the above inequality, we obtain

$$
\bar{x}^{\prime} P \bar{x}+\limsup _{k \rightarrow \infty, k \in K}\left(k x^{k^{\prime}} Q x^{k}\right) \leq 0
$$

Since, by the positive semidefiniteness of $Q, x^{k^{\prime}} Q x^{k} \geq 0$, we see that $\left\{x^{k^{\prime}} Q x^{k}\right\}_{k \in K}$ must converge to zero, for otherwise the left-hand side of the above inequality would be $+\infty$. Therefore, $\bar{x}^{\prime} Q \bar{x}=0$ and from our hypothesis we obtain $\bar{x}^{\prime} P \bar{x}>0$. This contradicts Eq. (2.5). Q.E.D.

Let us introduce now the augmented Lagrangian function

$$
L_{c}(x, \lambda)=f(x)+\lambda^{\prime} h(x)+\frac{c}{2}\|h(x)\|^{2}
$$

where $c$ is a scalar. This is the Lagrangian function for the problem

$$
\begin{aligned}
& \text { minimize } f(x)+\frac{c}{2}\|h(x)\|^{2} \\
& \text { subject to } h(x)=0
\end{aligned}
$$

which has the same local minima as our original problem of minimizing $f(x)$ subject to $h(x)=0$. The gradient and Hessian of $L_{c}$ with respect to $x$ are

$$
\begin{gathered}
\nabla_{x} L_{c}(x, \lambda)=\nabla f(x)+\nabla h(x)(\lambda+c h(x)) \\
\nabla_{\bar{x} x}^{2} L_{c}(x, \lambda)=\nabla^{2} f(x)+\sum_{i=1}^{m}\left(\lambda_{i}+c h_{i}(x)\right) \nabla^{2} h_{i}(x)+c \nabla h(x) \nabla h(x)^{\prime}
\end{gathered}
$$

In particular, if $x^{*}$ and $\lambda^{*}$ satisfy the sufficiency conditions of Prop. 3.2.1, we have

$$
\begin{aligned}
& \nabla_{x} L_{c}\left(x^{*}, \lambda^{*}\right)=\nabla f\left(x^{*}\right)+\nabla h\left(x^{*}\right)\left(\lambda^{*}+c h\left(x^{*}\right)\right)=\nabla_{x} L\left(x^{*}, \lambda^{*}\right)=0 \\
& \nabla_{\bar{x} x}^{2} L_{c}\left(x^{*}, \lambda^{*}\right)=\nabla^{2} f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla^{2} h_{i}\left(x^{*}\right)+c \nabla h\left(x^{*}\right) \nabla h\left(x^{*}\right)^{\prime} \\
& =\nabla_{\bar{x} x}^{2} L\left(x^{*}, \lambda^{*}\right)+c \nabla h\left(x^{*}\right) \nabla h\left(x^{*}\right)^{\prime}
\end{aligned}
$$

By the sufficiency condition (2.2), we have that $y^{\prime} \nabla_{\bar{x} x}^{2} L\left(x^{*}, \lambda^{*}\right) y>0$ for all $y \neq 0$ such that $y^{\prime} \nabla h\left(x^{*}\right) \nabla h\left(x^{*}\right)^{\prime} y=0$, so Lemma 3.2.1 implies that there exists a $\bar{c}$ such that

$$
\nabla_{\bar{x} x}^{2} L_{c}\left(x^{*}, \lambda^{*}\right): \text { positive definite, } \quad \forall c>\bar{c}
$$

---

![[ch3_and_4_p24_img6.jpeg]]

Figure 3.2.1. Illustration of how $x^{*}$ is an unconstrained minimum of the augmented Lagrangian $L_{c}\left(\cdot, \lambda^{*}\right)$ for sufficiently large $c$. Consider the two-dimensional problem

$$
\begin{aligned}
& \operatorname{minimize} f(x)=\frac{1}{2}\left(x_{1}^{2}-x_{2}^{2}\right)-x_{2} \\
& \text { subject to } x_{2}=0
\end{aligned}
$$

Here, $x^{*}=(0,0)$ is the unique global minimum, and $\lambda^{*}=1$ is a corresponding Lagrange multiplier. The augmented Lagrangian function is

$$
L_{c}\left(x, \lambda^{*}\right)=\frac{1}{2}\left(x_{1}^{2}-x_{2}^{2}\right)-x_{2}+\lambda^{*} x_{2}+\frac{c}{2} x_{2}^{2}=\frac{1}{2}\left(x_{1}^{2}-x_{2}^{2}\right)+\frac{c}{2} x_{2}^{2}
$$

and has $x^{*}=(0,0)$ as its unique unconstrained minimum for $c>1$. The figure shows the equal cost surfaces of $f$ (left side) and of $L_{c}\left(\cdot, \lambda^{*}\right)$ for $c=2$ (right side).

Using the sufficient optimality condition for unconstrained optimization (cf. Prop. 1.1.3), we conclude from Eqs. (2.6) and (2.7), that for $c>\bar{c}$. $x^{*}$ is an unconstrained local minimum of $L_{c}\left(\cdot, \lambda^{*}\right)$. In particular, there exist $\gamma>0$ and $\epsilon>0$ such that

$$
L_{c}\left(x, \lambda^{*}\right) \geq L_{c}\left(x^{*}, \lambda^{*}\right)+\frac{\gamma}{2}\left\|x-x^{*}\right\|^{2}, \quad \forall x \text { with }\left\|x-x^{*}\right\|<\epsilon
$$

Since for all $x$ with $h(x)=0$ we have $L_{c}\left(x, \lambda^{*}\right)=f(x)$, it follows that

$$
f(x) \geq f\left(x^{*}\right)+\frac{\gamma}{2}\left\|x-x^{*}\right\|^{2}, \quad \forall x \text { with } h(x)=0, \text { and }\left\|x-x^{*}\right\|<\epsilon
$$

---

Thus $x^{*}$ is a strict local minimum of $f$ over $h(x)=0$, providing an alternative proof of Prop. 3.2.1.

The preceding analysis also shows that we can try to minimize $f$ over $h(x)=0$ by computing an unconstrained minimum of the augmented Lagrangian $L_{c}\left(\cdot \lambda^{*}\right)$, as illustrated in Fig. 3.2.1. The difficulty here is that the Lagrange multiplier $\lambda^{*}$ is unknown, but it turns out that by using in place of $\lambda^{*}$, readily computable approximations to $\lambda^{*}$, we can obtain useful algorithms; see Section 4.2.

# 3.2.2 Sensitivity* 

Lagrange multipliers frequently have an interesting interpretation in specific practical contexts. In economic applications they can often be interpreted as prices, while in other problems they represent quantities with concrete physical meaning. It turns out that within our mathematical framework, they can be viewed as rates of change of the optimal cost as the level of constraint changes. This is fairly easy to show for the case of linear constraints as indicated in Fig. 3.2.2. The following proposition provides a formal proof for the case of nonlinear constraints.

Proposition 3.2.2: (Sensitivity Theorem) Let $f$ and $h_{i}$ be twice continuously differentiable and consider the family of problems

$$
\begin{aligned}
& \text { minimize } f(x) \\
& \text { subject to } h(x)=u
\end{aligned}
$$

parameterized by the vector $u \in \Re^{m}$. Assume that for $u=0$, the problem above has a local minimum $x^{*}$, which is regular and which together with its associated Lagrange multiplier vector $\lambda^{*}$ satisfies the second order sufficiency conditions of Prop. 3.2.1. Then there exists an open sphere $S$ centered at $u=0$ such that for every $u \in S$, there is an $x(u) \in \Re^{n}$ and a $\lambda(u) \in \Re^{m}$, which are a local minimum-Lagrange multiplier pair of problem (2.10). Furthermore, $x(\cdot)$ and $\lambda(\cdot)$ are continuously differentiable functions within $S$ and we have $x(0)=x^{*}$, $\lambda(0)=\lambda^{*}$. In addition, for all $u \in S$ we have

$$
\nabla p(u)=-\lambda(u)
$$

where $p(u)$ is the optimal cost parameterized by $u$, that is,

$$
p(u)=f(x(u))
$$

---

![[ch3_and_4_p26_img7.jpeg]]

Figure 3.2.2. Illustration of the sensitivity theorem for a problem involving a single linear constraint,

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } a^{\prime} x=b
\end{aligned}
$$

Here, $x^{*}$ is a local minimum and $\lambda^{*}$ is a corresponding Lagrange multiplier. If the level of constraint $b$ is changed to $b+\Delta b$, the minimum $x^{*}$ will change to $x^{*}+\Delta x$. Since $b+\Delta b=a^{\prime}\left(x^{*}+\Delta x\right)=a^{\prime} x^{*}+a^{\prime} \Delta x=b+a^{\prime} \Delta x$, we see that the variations $\Delta x$ and $\Delta b$ are related by

$$
a^{\prime} \Delta x=\Delta b
$$

Using the Lagrange multiplier condition $\nabla f\left(x^{*}\right)=-\lambda^{*} a$, the corresponding cost change can be written as

$$
\Delta \operatorname{cost}=f\left(x^{*}+\Delta x\right)-f\left(x^{*}\right)=\nabla f\left(x^{*}\right)^{\prime} \Delta x+o(\|\Delta x\|)=-\lambda^{*} a^{\prime} \Delta x+o(\|\Delta x\|)
$$

By combining the above two relations, we obtain $\Delta \operatorname{cost}=-\lambda^{*} \Delta b+o(\|\Delta x\|)$, so up to first order we have

$$
\lambda^{*}=-\frac{\Delta \cos t}{\Delta b}
$$

Thus, the Lagrange multiplier $\lambda^{*}$ gives the rate of optimal cost decrease as the level of constraint increases.

In the case where there are multiple constraints $a_{i}^{\prime} x=b_{i}, i=1, \ldots, m$, the preceding argument can be appropriately modified. In particular, we have

$$
\begin{aligned}
\Delta \text { cost } & =f\left(x^{*}+\Delta x\right)-f\left(x^{*}\right) \\
& =\nabla f\left(x^{*}\right)^{\prime} \Delta x+o(\|\Delta x\|) \\
& =-\sum_{i=1}^{m} \lambda_{i}^{*} a_{i}^{\prime} \Delta x+o(\|\Delta x\|)
\end{aligned}
$$

and $a_{i}^{\prime} \Delta x=\Delta b_{i}$ for all $i$, so we obtain $\Delta \operatorname{cost}=-\sum_{i=1}^{m} \lambda_{i}^{*} \Delta b_{i}+o(\|\Delta x\|)$.

---

Proof: Consider the system of equations

$$
\nabla f(x)+\nabla h(x) \lambda=0, \quad h(x)=u
$$

For each fixed $u$, this system represents $n+m$ equations with $n+m$ unknowns - the vectors $x$ and $\lambda$. For $u=0$ the system has the solution $\left(x^{*}, \lambda^{*}\right)$. The corresponding $(n+m) \times(n+m)$ Jacobian matrix with respect to $(x, \lambda)$ is given by

$$
J=\left(\begin{array}{cc}
\nabla^{2} f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla^{2} h_{i}\left(x^{*}\right) & \nabla h\left(x^{*}\right) \\
\nabla h\left(x^{*}\right)^{\prime} & 0
\end{array}\right)
$$

Let us show that $J$ is nonsingular. If it were not, some nonzero vector $\left(y^{\prime}, z^{\prime}\right)^{\prime}$ would belong to the nullspace of $J$, that is,

$$
\begin{gathered}
\left(\nabla^{2} f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla^{2} h_{i}\left(x^{*}\right)\right) y+\nabla h\left(x^{*}\right) z=0 \\
\nabla h\left(x^{*}\right)^{\prime} y=0
\end{gathered}
$$

Premultiplying Eq. (2.13) by $y^{\prime}$ and using Eq. (2.14), we obtain

$$
y^{\prime}\left(\nabla^{2} f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla^{2} h_{i}\left(x^{*}\right)\right) y=0
$$

In view of Eq. (2.14), it follows that $y=0$, for otherwise our second order sufficiency assumption would be violated [cf. Eq. (2.2)]. Since $y=0$, Eq. (2.13) yields $\nabla h\left(x^{*}\right) z=0$, which in view of the linear independence of the columns $\nabla h_{i}\left(x^{*}\right), i=1, \ldots, m$, of $\nabla h\left(x^{*}\right)$, yields $z=0 . . \mathrm{Thus}$, we obtain $y=0, z=0$, which is a contradiction. Hence, $J$ is nonsingular.

Returning now to the system (2.12), it follows from the nonsingularity of $J$ and the implicit function theorem (Prop. A. 25 in Appendix A) that for all $u$ in some open sphere $S$ centered at $u=0$, there exist $x(u)$ and $\lambda(u)$ such that $x(0)=x^{*}, \lambda(0)=\lambda^{*}$, the functions $x(\cdot)$ and $\lambda(\cdot)$ are continuously differentiable, and

$$
\begin{gathered}
\nabla f(x(u))+\nabla h(x(u)) \lambda(u)=0 \\
h(x(u))=u
\end{gathered}
$$

For $u$ sufficiently close to $u=0$, the vectors $x(u)$ and $\lambda(u)$ satisfy the second order sufficiency conditions for problem (2.10), since they satisfy them by assumption for $u=0$. This is straightforward to verify by using our continuity assumptions. [If it were not true, there would exist a sequence $\left\{u^{k}\right\}$

---

![[ch3_and_4_p28_img8.jpeg]]

Figure 3.2.3. Illustration of the primal function $p(u)=f(x(u))$ for the twodimensional problem of Fig. 3.2.1

$$
\begin{aligned}
& \operatorname{minimize} f(x)=\frac{1}{2}\left(x_{1}^{2}-x_{2}^{2}\right)-x_{2} \\
& \text { subject to } h(x)=x_{2}=0
\end{aligned}
$$

Here,

$$
p(u)=\min _{h(x)=u} f(x)=-\frac{1}{2} u^{2}-u
$$

and we have $\lambda^{*}=-\nabla p(0)=1$ consistently with the sensitivity theorem.
with $u^{k} \rightarrow 0$, and a sequence $\left\{y^{k}\right\}$ with $\left\|y^{k}\right\|=1$ and $\nabla h\left(x\left(u^{k}\right)\right)^{\prime} y^{k}=0$ for all $k$, such that

$$
y^{k^{\prime}} \nabla_{x x}^{2} L\left(x\left(u^{k}\right), \lambda\left(u^{k}\right)\right) y^{k} \leq 0, \quad \forall k
$$

By taking the limit along a convergent subsequence of $\left\{y^{k}\right\}$, we would obtain a contradiction of the second order sufficiency condition at $\left(x^{*}, \lambda^{*}\right)$; cf. Eq. (2.2).] Hence, $x(u)$ and $\lambda(u)$ are a local minimum-Lagrange multiplier pair for problem (2.10).

There remains to show that $\nabla p(u)=\nabla_{u}\{f(x(u))\}=-\lambda(u)$. By multiplying Eq. (2.15) by $\nabla x(u)$, we obtain

$$
\nabla x(u) \nabla f(x(u))+\nabla x(u) \nabla h(x(u)) \lambda(u)=0
$$

By differentiating the relation $h(x(u))=u$, it follows that

$$
I=\nabla_{u}\{h(x(u))\}=\nabla x(u) \nabla h(x(u))
$$

where $I$ is the $m \times m$ identity matrix. Finally, by using the chain rule, we have

$$
\nabla p(u)=\nabla_{u}\{f(x(u))\}=\nabla x(u) \nabla f(x(u))
$$

---

Combining the above three relations, we obtain

$$
\nabla p(u)+\lambda(u)=0
$$

and the proof is complete. Q.E.D.
The function $p(u)=f(x(u))$ in the sensitivity theorem is called the primal function and plays an important in contexts other than sensitivity. When we discuss duality theory later, we will see that it is closely related to the, so-called, dual function. Within that context, the relation $\nabla p(0)=$ $-\lambda^{*}$, illustrated in Fig. 3.2.3, will also turn out to be significant.

Note that the regularity of $x^{*}$ is essential for the primal function to be defined within a sphere centered at 0 . For example, for the problem of minimizing $\frac{1}{2}\left(x_{1}^{2}-x_{2}^{2}\right)-x_{2}$ subject to $\left(x_{2}\right)^{2}=0$, which is equivalent to the problem of Fig. 3.2.3, the primal function is undefined for $u<0$.

# E XERCISES 

2.1

Consider the problem

$$
\begin{aligned}
& \text { minimize } x_{1}+x_{2} \\
& \text { subject to } x_{1}^{2}+x_{2}^{2}=2
\end{aligned}
$$

(a) Use the second order sufficiency conditions to verify that it has a unique global minimum and a unique global maximum.
(b) Calculate the primal function corresponding to the minimum and verify that its gradient is related to the Lagrange multiplier as specified by the sensitivity theorem.

## 2.2

Consider the problem of Example 1.3 in Section 3.1. Calculate the primal function corresponding to the minimum and verify that its gradient is related to the Lagrange multiplier as specified by the sensitivity theorem.

---

# 2.3 (Hessian of the Primal Function) 

Under the assumptions of the sensitivity theorem, show that for every scalar $c$ for which the matrix

$$
A_{c}(u)=\nabla_{x x}^{2} L(x(u), \lambda(u))+c \nabla h(x(u)) \nabla h(x(u))^{\prime}
$$

is invertible, we have

$$
\nabla^{2} p(u)=\nabla h(x(u))^{\prime} A_{c}(u)^{-1} \nabla h(x(u))-c I
$$

Hint: Differentiate Eqs. (2.15) and (2.18), and use Eq. (2.17).

### 3.3 INEQUALITY CONSTRAINTS

Consider a problem involving both equality and inequality constraints

$$
\begin{array}{cl}
\operatorname{minimize} & f(x) \\
\text { subject to } & h_{1}(x)=0, \ldots, h_{m}(0)=0 \\
& g_{1}(x) \leq 0, \ldots, g_{r}(x) \leq 0
\end{array}
$$

where $f, h_{i}, g_{j}$ are continuously differentiable functions from $\Re^{n}$ to $\Re$. Using the notation

$$
h=\left(h_{1}, \ldots, h_{m}\right), \quad g=\left(g_{1}, \ldots, g_{r}\right)
$$

the problem can be written as

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } h(x)=0, \quad g(x) \leq 0
\end{aligned}
$$

For any feasible point $x$, the set of active inequality constraints is denoted by

$$
A(x)=\left\{j \mid g_{j}(x)=0\right\}
$$

If $j \notin A(x)$, we say that the $j$ th constraint is inactive at $x$. We note that if $x^{*}$ is a local minimum of problem (ICP), then $x^{*}$ is also a local minimum for a problem identical to (ICP) except that the inactive constraints at $x^{*}$ have been discarded. Thus, in effect, inactive constraints at $x^{*}$ don't matter; they can be ignored in the statement of optimality conditions.

On the other hand, at a local minimum, active inequality constraints can be treated to a large extent as equalities. In particular, if $x^{*}$ is a local

---

minimum of the inequality constrained problem (ICP), then $x^{*}$ is also a local minimum for the equality constrained problem
minimize $f(x)$
subject to $h_{1}(x)=0, \ldots, h_{m}(x)=0, \quad g_{j}(x)=0, \quad \forall j \in A\left(x^{*}\right)$.
Thus, if $x^{*}$ is regular for the latter problem, there exist Lagrange multipliers $\lambda_{1}^{*}, \ldots, \lambda_{m}^{*}$, and $\mu_{j}^{*}, j \in A\left(x^{*}\right)$ such that

$$
\nabla f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla h_{i}\left(x^{*}\right)+\sum_{j \in A\left(x^{*}\right)} \mu_{j}^{*} \nabla g_{j}\left(x^{*}\right)=0
$$

Assigning zero Lagrange multipliers to the inactive constraints, we obtain

$$
\begin{gathered}
\nabla f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla h_{i}\left(x^{*}\right)+\sum_{j=1}^{r} \mu_{j}^{*} \nabla g_{j}\left(x^{*}\right)=0 \\
\mu_{j}^{*}=0, \quad \forall j \notin A\left(x^{*}\right)
\end{gathered}
$$

which may be viewed as an analog of the first order optimality condition for the equality constrained problem.

There is one more important fact about the Lagrange multipliers $\mu_{j}^{*}$; they are nonnegative. For a geometric illustration, see Fig. 3.3.1. For an algebraic argument, note that if the $j$ th constraint is relaxed to $g_{j}(x) \leq u_{j}$, where $u_{j}>0$, the optimal cost will tend to decrease because the constraint set will become larger. Therefore, the sensitivity interpretation of the Lagrange multiplier $\left[\mu_{j}^{*}=-\left(\Delta\right.\right.$ cost due to $\left.\left.u_{j}\right) / u_{j}\right.$, cf. Prop. 3.2.2] suggests that $\mu_{j}^{*} \geq 0$.
![[ch3_and_4_p31_img9.jpeg]]

Figure 3.3.1 Illustration of the nonnegativity of the Lagrange multiplier for a problem with a single inequality constraint. If the constraint is inactive, then $\mu^{*}=0$. Otherwise, $\nabla f\left(x^{*}\right)$ is normal to the constraint surface and points to the inside of the constraint set, while $\nabla g\left(x^{*}\right)$ is normal to the constraint surface and points to the outside of the constraint set. Thus $\nabla f\left(x^{*}\right)$ and $\nabla g\left(x^{*}\right)$ are collinear and have opposite signs, implying that the Lagrange multiplier is nonnegative.

---

# 3.3.1 Necessary Conditions 

The preceding informal discussion will now be made rigorous. We first extend the notion of regularity that we used for equality constrained problems. A feasible vector $x$ is said to be regular if the equality constraint gradients $\nabla h_{i}(x), i=1, \ldots, m$, and the active inequality constraint gradients $\nabla g_{j}(x), j \in A(x)$, are linearly independent.

The following proposition generalizes the Lagrange multiplier theorem of Section 3.1. It will be stated in terms of the Lagrangian function defined by

$$
L(x, \lambda, \mu)=f(x)+\sum_{i=1}^{m} \lambda_{i} h_{i}(x)+\sum_{j=1}^{r} \mu_{j} g_{j}(x)
$$

Proposition 3.3.1: (Kuhn-Tucker Necessary Conditions) Let $x^{*}$ be a local minimum of $f$ subject to $h(x)=0, g(x) \leq 0$, and assume that $x^{*}$ is regular. Then there exist unique Lagrange multiplier vectors $\lambda^{*}=\left(\lambda_{1}^{*}, \ldots, \lambda_{m}^{*}\right), \mu^{*}=\left(\mu_{1}^{*}, \ldots, \mu_{r}^{*}\right)$, such that

$$
\begin{gathered}
\nabla_{x} L\left(x^{*}, \lambda^{*}, \mu^{*}\right)=0 \\
\mu_{j}^{*} \geq 0, \quad j=1, \ldots, r \\
\mu_{j}^{*}=0, \quad \forall j \notin A\left(x^{*}\right)
\end{gathered}
$$

where $A\left(x^{*}\right)$ is the set of active constraints at $x^{*}$ [cf. Eq. (3.2)]. If in addition $f, h$, and $g$ are twice continuously differentiable, there holds

$$
y^{\prime} \nabla_{x x}^{2} L\left(x^{*}, \lambda^{*}, \mu^{*}\right) y \geq 0, \quad \text { for all } y \in V\left(x^{*}\right)
$$

where

$$
V\left(x^{*}\right)=\left\{y \mid \nabla h_{i}\left(x^{*}\right)^{\prime} y=0, i=1, \ldots, m, \nabla g_{j}\left(x^{*}\right)^{\prime} y=0, j \in A\left(x^{*}\right)\right\}
$$

Proof: All the assertions of the proposition follow from the preceding discussion, and the first and second order optimality conditions for the equality constrained problem (3.3) (cf. Prop. 3.1.1), except for the assertion $\mu_{j}^{*} \geq 0$ for $j \in A\left(x^{*}\right)$. We give a proof of this assertion using the penalty approach of Section 3.1. We introduce the functions

$$
g_{j}^{+}(x)=\max \left\{0, g_{j}(x)\right\}, \quad j=1, \ldots, r
$$

---

and for each $k=1,2, \ldots$, the "penalized" problem

$$
\begin{aligned}
& \operatorname{minimize} F^{k}(x) \equiv f(x)+\frac{k}{2}\|h(x)\|^{2}+\frac{k}{2} \sum_{j=1}^{r}\left(g_{j}^{+}(x)\right)^{2}+\frac{1}{2}\left\|x-x^{*}\right\|^{2} \\
& \text { subject to } x \in S \text {, } \\
& \text { where } S=\left\{x \mid\left\|x-x^{*}\right\| \leq \epsilon\right\}, \text { and } \epsilon>0 \text { is such that } f\left(x^{*}\right) \leq f(x) \text { for } \\
& \text { all feasible } x \text { with } x \in S \text {. Note that the function }\left(g_{j}^{+}(x)\right)^{2} \text { is continuously } \\
& \text { differentiable with gradient } 2 g_{j}^{+}(x) \nabla g_{j}(x) \text {. If } x^{k} \text { minimizes } F^{k}(x) \text { over } S \text {, a } \\
& \text { similar argument to the one used for the equality constraint case in Section } \\
& 3.1 \text { shows that } x^{k} \rightarrow x^{*}, \text { and that the Lagrange multipliers } \lambda_{i}^{*} \text { and } \mu_{j}^{*} \text { are } \\
& \text { given by } \\
& \lambda_{i}^{*}=\lim _{k \rightarrow \infty} k h_{i}\left(x^{k}\right), \quad i=1, \ldots, m, \\
& \mu_{j}^{*}=\lim _{k \rightarrow \infty} k g_{j}^{+}\left(x^{k}\right), \quad j=1, \ldots, r .
\end{aligned}
$$

Since $g_{j}^{+}\left(x^{k}\right) \geq 0$, we obtain $\mu_{j}^{*} \geq 0$ for all $j$. Q.E.D.

# Example 3.1 

Consider the problem

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2}\left(x_{1}^{2}+x_{2}^{2}+x_{3}^{2}\right) \\
& \text { subject to } x_{1}+x_{2}+x_{3} \leq-3
\end{aligned}
$$

Then for a local minimum $x^{*}$, the first order necessary condition [cf. Eq. (3.4)] yields

$$
\begin{aligned}
& x_{1}^{*}+\mu^{*}=0 \\
& x_{2}^{*}+\mu^{*}=0 \\
& x_{3}^{*}+\mu^{*}=0
\end{aligned}
$$

There are two possibilities:
(a) The constraint is inactive,

$$
x_{1}^{*}+x_{2}^{*}+x_{3}^{*}<-3
$$

in which case $\mu^{*}=0$. Then, we obtain $x_{1}^{*}=x_{2}^{*}=x_{3}^{*}=0$, which contradicts the inequality $x_{1}^{*}+x_{2}^{*}+x_{3}^{*}<-3$. Hence, this possibility is excluded.
(b) The constraint is active,

$$
x_{1}^{*}+x_{2}^{*}+x_{3}^{*}=-3
$$

Then we obtain $x_{1}^{*}=x_{2}^{*}=x_{3}^{*}=-1$ and $\mu^{*}=1$, which satisfy all the necessary conditions for a local minimum. Since every point is regular it follows that $x^{*}=(-1,-1,-1)$ is the unique candidate for a local minimum. We may proceed further and verify that $x^{*}$ and $\mu^{*}$ satisfy the second order sufficiency conditions for optimality (to be provided shortly) and conclude that $x^{*}$ is indeed the unique local minimum for this problem.

---

# 3.3.2 Proof Based on Conversion to the Equality Case* 

We now provide an alternative proof of the necessary conditions of Prop. 3.3.1, based on converting the inequality constrained problem (ICP) into an equality constrained problem and then using the corresponding necessary conditions of the previous section. This approach is straightforward and is also useful in other contexts. It is not as general, however, as our earlier approach in this section, because it requires twice differentiability of the problem functions.

Consider the equality constrained problem

$$
\begin{aligned}
& \text { minimize } f(x) \\
& \text { subject to } \quad h_{1}(x)=0, \ldots, h_{m}(x)=0 \\
& g_{1}(x)+z_{1}^{2}=0, \ldots, g_{r}(x)+z_{r}^{2}=0
\end{aligned}
$$

where we have introduced additional variables $z_{1}, \ldots, z_{r}$.
Let $x^{*}$ be a local minimum for our original problem (ICP). Then $\left(x^{*}, z^{*}\right)$ is a local minimum for problem (3.10), where $z^{*}=\left(z_{1}^{*}, \ldots, z_{r}^{*}\right)$,

$$
z_{j}^{*}=\left(-g_{j}\left(x^{*}\right)\right)^{1 / 2}, \quad j=1, \ldots, r
$$

It is straightforward to verify that $x^{*}$ is regular for problem (ICP), if and only if ( $x^{*}, z^{*}$ ) is regular for problem (3.10); we leave this as an exercise for the reader. By applying the first order necessary conditions for equality constraints (Prop. 3.1.1) to problem (3.10), we see that there exist Lagrange multipliers $\lambda_{1}^{*}, \ldots, \lambda_{m}^{*}, \mu_{1}^{*}, \ldots, \mu_{r}^{*}$ such that

$$
\begin{gathered}
\nabla f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} \nabla h_{i}\left(x^{*}\right)+\sum_{j=1}^{r} \mu_{j}^{*} \nabla g_{j}\left(x^{*}\right)=0 \\
2 \mu_{j}^{*} z_{j}^{*}=0, \quad j=1, \ldots, r
\end{gathered}
$$

Since $z_{j}^{*}=\left(-g_{j}\left(x^{*}\right)\right)^{1 / 2}>0$ for $j \notin A\left(x^{*}\right)$, the last equation can also be written as

$$
\mu_{j}^{*}=0, \quad \forall j \notin A\left(x^{*}\right)
$$

Thus, to prove Prop. 3.3.1, there remains to show the nonnegativity condition $\mu_{j}^{*} \geq 0$ and the second order condition (3.8). For this purpose we use the second order necessary condition for the equivalent equality constrained problem (3.10). It yields [cf. Eq. (1.24) in Section 3.1]
$\left(y^{\prime} \quad w^{\prime}\right)\left(\begin{array}{ccccc}\nabla_{x x}^{2} L\left(x^{*}, \lambda^{*}, \mu^{*}\right) & & 0 & & \\ & 2 \mu_{1}^{*} & 0 & \ldots & 0 \\ 0 & 0 & 2 \mu_{2}^{*} & \ldots & 0 \\ & \vdots & \vdots & \vdots & \vdots \\ 0 & 0 & \ldots & 2 \mu_{r}^{*}\end{array}\right)\binom{y}{w} \geq 0$,

---

for all $y \in \Re^{n}, w \in \Re^{r}$ satisfying

$$
\nabla h\left(x^{*}\right)^{\prime} y=0, \quad \nabla g_{j}\left(x^{*}\right)^{\prime} y+2 z_{j}^{*} w_{j}=0, \quad j=1, \ldots, r
$$

We will now try different pairs $(y, w)$ satisfying Eq. (3.15) and obtain the desired results from Eq. (3.14). First let us select $w$ by

$$
w_{j}= \begin{cases}0 & \text { if } j \in A\left(x^{*}\right) \\ -\frac{\nabla g_{j}\left(x^{*}\right)^{\prime} y}{2 z_{j}^{*}} & \text { if } j \notin A\left(x^{*}\right)\end{cases}
$$

Then $\nabla g_{j}\left(x^{*}\right)^{\prime} y+2 z_{j}^{*} w_{j}=0$ for all $j \notin A\left(x^{*}\right), z_{j}^{*} w_{j}=0$ for all $j \in A\left(x^{*}\right)$, and by taking into account the fact $\mu_{j}^{*}=0$ for $j \notin A\left(x^{*}\right)$ [cf. Eq. (3.13)], we also have $\mu_{j}^{*} w_{j}=0$ for all $j$. Thus we obtain from Eqs. (3.14) and (3.15)

$$
y^{\prime} \nabla_{z x}^{2} L\left(x^{*}, \lambda^{*}, \mu^{*}\right) y \geq 0
$$

for all $y$ such that $\nabla h\left(x^{*}\right)^{\prime} y=0$ and $\nabla g_{j}\left(x^{*}\right)^{\prime} y=0$ for all $j \in A\left(x^{*}\right)$, thereby verifying the second order condition (3.8).

Next let us select, for every $j \in A\left(x^{*}\right)$, a vector $(y, w)$ with $y=0$, $w_{j} \neq 0, w_{k}=0$ for all $k \neq j$. Such a vector satisfies the condition of Eq. (3.15). By using such a vector in Eq. (3.14), we obtain $2 \mu_{j}^{*} w_{j}^{2} \geq 0$, and

$$
\mu_{j}^{*} \geq 0, \quad \forall j \in A\left(x^{*}\right)
$$

Thus, the conversion to the equivalent equality constrained problem (3.10) has yielded all the Kuhn-Tucker conditions for the inequality constrained problem.

# 3.3.3 Second Order Sufficiency Conditions and Sensitivity* 

The transformation to an equality constrained problem described above may also be used to derive sufficiency conditions for the inequality constrained problem (ICP).

Proposition 3.3.2: (Second Order Sufficiency Conditions) Assume that $f, h$, and $g$ are twice continuously differentiable, and let $x^{*} \in \Re^{n}, \lambda^{*} \in \Re^{m}$, and $\mu^{*} \in \Re^{r}$ satisfy

$$
\nabla_{x} L\left(x^{*}, \lambda^{*}, \mu^{*}\right)=0, \quad h\left(x^{*}\right)=0, \quad g\left(x^{*}\right) \leq 0
$$

---

$$
\begin{array}{ll}
\mu_{j}^{*} \geq 0, & j=1, \ldots, r \\
\mu_{j}^{*}=0, & \forall j \notin A\left(x^{*}\right) \\
y^{\prime} \nabla_{x x}^{2} L\left(x^{*}, \lambda^{*}, \mu^{*}\right) y>0, & \text { for all } y \neq 0 \text { with } y \in V\left(x^{*}\right)
\end{array}
$$

where

$$
V\left(x^{*}\right)=\left\{y \mid \nabla h_{i}\left(x^{*}\right)^{\prime} y=0, i=1, \ldots, m, \nabla g_{j}\left(x^{*}\right)^{\prime} y=0, j \in A\left(x^{*}\right)\right\}
$$

Assume also that

$$
\mu_{j}^{*}>0, \quad \forall j \in A\left(x^{*}\right)
$$

Then $x^{*}$ is a strict local minimum of $f$ subject to $h(x)=0, g(x) \leq 0$.

Proof: (Abbreviated) It is straightforward to verify using our assumptions, that for the equivalent equality constrained problem (3.10), the second order sufficiency conditions are satisfied by $\left(x^{*}, z_{1}^{*}, \ldots, z_{r}^{*}\right), \lambda^{*}, \mu^{*}$, where $z_{j}^{*}=\left(-g_{j}\left(x^{*}\right)\right)^{1 / 2}$. Hence $\left(x^{*}, z_{1}^{*}, \ldots, z_{r}^{*}\right)$ is a strict local minimum for problem (3.10) and the result follows. Q.E.D.

The condition (3.21) is known as the strict complementarity assumption. Exercise 3.7 provides a sharper version of the preceding sufficiency theorem, where the strict complementarity assumption is somewhat weakened. However, this assumption cannot be entirely discarded as the following example shows.

# Example 3.2 

Consider the two-dimensional problem

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2}\left(x_{1}^{2}-x_{2}^{2}\right) \\
& \text { subject to } x_{2} \leq 0
\end{aligned}
$$

It can be verified that $x^{*}=(0,0)$ and $\mu^{*}=0$ satisfy all the conditions of Prop. 3.3.2 except for the strict complementarity assumption. It is seen that $x^{*}$ is not a local minimum, because the cost can be decreased by moving from $x^{*}$ along the feasible direction $(0,-1)$. This direction is not in the subspace $V\left(x^{*}\right)$ and is not accounted for in the condition (3.19). The problem here is that the constraint is active but degenerate, in the sense that its Lagrange multiplier is zero (see Exercise 3.7).

The transformation to an equality constrained problem together with the corresponding sensitivity result of Prop. 3.2.2 also yields the following proposition. The proof is left for the reader.

---

Proposition 3.3.3: (Sensitivity) Let $f, h$, and $g$ be twice continuously differentiable and consider the family of problems

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } h(x)=u, \quad g(x) \leq v
\end{aligned}
$$

parameterized by the vectors $u \in \Re^{m}$ and $v \in \Re^{r}$. Assume that for $(u, v)=(0,0)$ this problem has a local minimum $x^{*}$, which is regular and which together with its associated Lagrange multiplier vectors $\lambda^{*}$ and $\mu^{*}$, satisfies the second order sufficiency conditions of Prop. 3.3.2. Then there exists an open sphere $S$ centered at $(u, v)=(0,0)$ such that for every $(u, v) \in S$ there is an $x(u, v) \in \Re^{n}$ and $\lambda(u, v) \in \Re^{m}$, $\mu(u, v) \in \Re^{r}$, which are a local minimum and associated Lagrange multiplier vectors of problem (3.22). Furthermore, $x(\cdot), \lambda(\cdot)$, and $\mu(\cdot)$ are continuously differentiable in $S$ and we have $x(0,0)=x^{*}, \lambda(0,0)=$ $\lambda^{*}, \mu(0,0)=\mu^{*}$. In addition, for all $(u, v) \in S$, there holds

$$
\begin{aligned}
& \nabla_{u} p(u, v)=-\lambda(u, v) \\
& \nabla_{v} p(u, v)=-\mu(u, v)
\end{aligned}
$$

where $p(u, v)$ is the optimal cost parameterized by $(u, v)$,

$$
p(u, v)=f(x(u, v))
$$

# E XERCISES 

## 3.1

Solve the two-dimensional problem

$$
\begin{aligned}
& \operatorname{minimize}(x-a)^{2}+(y-b)^{2}+x y \\
& \text { subject to } 0 \leq x \leq 1, \quad 0 \leq y \leq 1
\end{aligned}
$$

for all possible values of the scalars $a$ and $b$.

---

Given a vector $y$, consider the problem

$$
\begin{aligned}
& \text { maximize } y^{\prime} x \\
& \text { subject to } x^{\prime} Q x \leq 1
\end{aligned}
$$

where $Q$ is a positive definite symmetric matrix. Show that the optimal value is $\sqrt{y^{\prime} Q^{-1} y}$ and use this fact to establish the inequality

$$
\left(x^{\prime} y\right)^{2} \leq\left(x^{\prime} Q x\right)\left(y^{\prime} Q^{-1} y\right)
$$

# 3.3 

Solve Exercise 1.14 of Section 2.1 using Lagrange multipliers.

## 3.4

Solve Exercise 1.15 of Section 2.1 using Lagrange multipliers.

### 3.5 (Constraint Qualification for Inequality Constraints)

Consider the problem

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } g(x) \leq 0
\end{aligned}
$$

For a feasible point $x$, let $D(x)$ denote the set of all feasible directions at $x$ defined by

$$
D(x)=\{d \mid d \neq 0, \text { and for some } \bar{\alpha}>0, g(x+\alpha d) \leq 0 \text { for all } \alpha \in[0, \bar{\alpha}]\}
$$

and denote by $\bar{D}(x)$ the closure of $D(x)$. Let $x^{*}$ be a local minimum for problem (3.23). Show that:
(a)

$$
\nabla f\left(x^{*}\right)^{\prime} d \geq 0, \quad \forall d \in \bar{D}\left(x^{*}\right)
$$

(b) If we have

$$
\bar{D}\left(x^{*}\right)=\left\{d \mid \nabla g_{j}\left(x^{*}\right)^{\prime} d \leq 0, \forall j \in A\left(x^{*}\right)\right\}
$$

then there exists a Lagrange multiplier vector $\mu^{*} \geq 0$ such that

$$
\nabla f\left(x^{*}\right)+\nabla g\left(x^{*}\right) \mu^{*}=0
$$

---

$$
\mu_{j}^{*}=0, \quad \forall j \notin A\left(x^{*}\right)
$$

[Condition (3.24) is called a constraint qualification.] Hint: Use part (a) and Farkas' Lemma [Prop. B.16(d) in Appendix B].
(c) The constraint qualification holds if any one of the following conditions holds:

1. The functions $g_{j}$ are linear.
2. [MaF67] There exists a vector $d$ such that

$$
\nabla g_{j}\left(x^{*}\right)^{\prime} d<0, \quad \forall j \in A\left(x^{*}\right)
$$

3. [Sla50] The functions $g_{j}$ are convex and there exists a vector $\bar{x}$ such that

$$
g_{j}(\bar{x})<0, \quad \forall j \in A\left(x^{*}\right)
$$

Hint: Let $d=\bar{x}-x^{*}$ and use condition 2.
4. The gradients $\nabla g_{j}\left(x^{*}\right), j \in A\left(x^{*}\right)$, are linearly independent. Hint: Let $d$ be such that $\nabla g_{j}\left(x^{*}\right)^{\prime} d=-1$ for all $j \in A\left(x^{*}\right)$ and use condition 2 .

# 3.6 (Minimax Problems) 

Consider the problem

$$
\begin{aligned}
& \operatorname{minimize} \max \left\{g_{1}(x), \ldots, g_{r}(x)\right\} \\
& \text { subject to } x \in \Re^{n}
\end{aligned}
$$

where $g_{j}: \Re^{n} \mapsto \Re$ are continuously differentiable. Show that if $x^{*}$ is a local minimum, there exists a vector $\mu^{*}=\left(\mu_{1}^{*}, \ldots, \mu_{r}^{*}\right)$ such that

$$
\begin{gathered}
\sum_{j=1}^{r} \mu_{j}^{*} \nabla g_{j}\left(x^{*}\right)=0 . \quad \mu^{*} \geq 0, \quad \sum_{j=1}^{r} \mu_{j}^{*}=1 \\
\mu_{j}^{*}=0, \quad \text { if } \quad g_{j}\left(x^{*}\right)<\max \left\{g_{1}\left(x^{*}\right), \ldots, g_{r}\left(x^{*}\right)\right\}
\end{gathered}
$$

Hint: Consider the equivalent problem $\min \left\{z \mid g_{j}(x) \leq z, j=1, \ldots, r\right\}$. Use also the results of Exercise 3.5(c).

## 3.7 (Sufficiency Without Strict Complementarity)

Show that the sufficiency result of Prop. 3.3.2 holds if the strict complementarity assumption (3.21) is removed, but the assumption on the Hessian of the Lagrangian is strengthened to read

$$
y^{\prime} \nabla_{x x}^{2} L\left(x^{*}, \lambda^{*}, \mu^{*}\right) y>0
$$

for all $y \neq 0$ satisfying

$$
\begin{aligned}
\nabla h_{i}\left(x^{*}\right)^{\prime} y=0, \quad i=1, \ldots, m \\
\nabla g_{j}\left(x^{*}\right)^{\prime} y=0, \quad \forall j \in A\left(x^{*}\right) \text { with } \mu_{j}^{*}>0 \\
\nabla g_{j}\left(x^{*}\right)^{\prime} y \leq 0, \quad \forall j \in A\left(x^{*}\right) \text { with } \mu_{j}^{*}=0
\end{aligned}
$$

Hint: Extend the proof of Prop. 3.2.1.

---

# 3.4 LINEAR CONSTRAINTS* 

Problems with linear constraints have a remarkable property. They possess Lagrange multipliers always, even for local minima which are not regular. This is a fairly deep consequence of the geometry of polyhedral sets. A key result is Farkas' lemma, given as Prop. B.16(d) in Appendix B. It says that the cone $C$ generated by $r$ vectors $a_{1}, \ldots, a_{r}$,

$$
C=\left\{x \mid x=\sum_{j=1}^{r} \mu_{j} a_{j}, \mu_{j} \geq 0, j=1, \ldots, r\right\}
$$

is closed and consists of all vectors $x$ satisfying

$$
x^{\prime} y \leq 0, \quad \text { for all } y \text { such that } a_{j}^{\prime} y \leq 0, \forall j=1, \ldots, r
$$

The lemma is illustrated in Fig. 3.4.1. An examination of its proof, given in Appendix B and outlined in Fig. 3.4.1, shows that once the cone $C$ is shown to be closed, the rest is relatively easy; however, showing that $C$ is closed is nontrivial (it is not generally true that the linear transformation of a closed convex set or the vector sum of a finite number of closed convex sets is a closed convex set). Unfortunately, many optimization books treat the closure property of $C$ as obvious, thus diverting attention from the mathematically deepest aspect of Farkas' lemma and the subsequent Lagrange multiplier theorem.

Consider the problem

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } a_{i}^{\prime} x \leq b_{j}, \quad j=1, \ldots, r
\end{aligned}
$$

where $a_{j}$ and $b_{j}$ are given vectors and scalars, respectively, and $f: \mathbb{R}^{n} \rightarrow \mathbb{R}$ is continuously differentiable. The following necessary condition for optimality is a straightforward consequence of Farkas' lemma.

Proposition 3.4.1: (Lagrange Multiplier Theorem - Linear Constraints) Let $x^{*}$ be a local minimum of problem (4.3). Then there exist scalars $\mu_{1}^{*}, \ldots, \mu_{r}^{*}$ with $\mu_{j}^{*} \geq 0, j=1, \ldots, r$, such that

$$
\begin{gathered}
\nabla f\left(x^{*}\right)+\sum_{j=1}^{r} \mu_{j}^{*} a_{j}=0 \\
\mu_{j}^{*}=0, \quad \forall j \notin A\left(x^{*}\right)
\end{gathered}
$$

where $A\left(x^{*}\right)$ is the set $\left\{j \mid a^{\prime} x^{*}=b_{j}\right\}$ of active constraints at $x^{*}$.

---

![[ch3_and_4_p41_img10.jpeg]]

Figure 3.4.1. Geometrie interpretation and outline of the proof of Farkas' lemma. Consider the cone $C$ generated by $r$ vectors $a_{1}, \ldots, a_{r}$,

$$
C=\left\{x \mid x=\sum_{j=1}^{r} \mu_{j} a_{j}, \mu_{j} \geq 0, j=1, \ldots, r\right\}
$$

and the cone

$$
C^{\perp}=\left\{y \mid a_{j}^{\prime} y \leq 0, j=1, \ldots, r\right\}
$$

Farkas' lemma says that if $x$ is such that

$$
x^{\prime} y \leq 0, \quad \forall y \in C^{\perp}
$$

then $x$ must belong to $C$ (the reverse is also clearly true). The proof consists of showing that C is closed this is the hard part shown in Prop. B.16(b) of Appendix B] and then considering the projection $\hat{x}$ on $C$ of a vector $x$ satisfying Eq. (*). From the geometry of projections we have

$$
x^{\prime}(x-\hat{x})=\|x-\hat{x}\|^{2}
$$

while

$$
(x-\hat{x})^{\prime} a_{j} \leq 0, \quad j=1, \ldots, r
$$

Hence, $(x-\hat{x}) \in C^{\perp}$, and using the hypothesis (*), we have

$$
x^{\prime}(x-\hat{x}) \leq 0
$$

Combining Eqs. ( $* *$ ) and ( $* * *$ ), we obtain $x=\hat{x}$, so $x \in C$.

Proof: See Fig. 3.4.2.
We finally note that linear equality constraints of the form $e_{i}^{\prime} x=$ $d_{i}, i=1, \ldots, m$, can be handled by conversion into the two inequality constraints $e_{i}^{\prime} x-d_{i} \leq 0$ and $-e_{i}^{\prime} x+d_{i} \leq 0$. Application of Prop. 3.4.1

---

![[ch3_and_4_p42_img11.jpeg]]

Figure 3.4.2. Proof of Prop. 3.4.1. We note that the local minimum $x^{*}$ is also a local minimum for the problem

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } a_{j}^{\prime} x \leq b_{j}, \quad j \in A\left(x^{*}\right)
\end{aligned}
$$

where $A\left(x^{*}\right)$ is the set of active constraints at $x^{*}$. By Prop. 2.1.2 in Section 2.1, we have that

$$
\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right) \geq 0, \quad \forall x \text { such that } a_{j}^{\prime} x \leq b_{j}, j \in A\left(x^{*}\right)
$$

Since a constraint $a_{j}^{\prime} x \leq b_{j}, j \in A\left(x^{*}\right)$ can also be expressed as $a_{j}^{\prime}\left(x-x^{*}\right) \leq 0$, the preceding condition is equivalent to

$$
\nabla f\left(x^{*}\right)^{\prime} y \geq 0, \quad \forall y \text { such that } a_{j}^{\prime} y \leq 0, j \in A\left(x^{*}\right)
$$

From Farkas' lemma, it follows that $-\nabla f\left(x^{*}\right)$ is in the cone generated by the vectors $a_{j}, j \in A\left(x^{*}\right)$ or equivalently, that there exist nonnegative scalars $\mu_{j}^{*}$, $j \in A\left(x^{*}\right)$, such that

$$
\nabla f\left(x^{*}\right)+\sum_{j \in A\left(x^{*}\right)} \mu_{j}^{*} a_{j}=0
$$

By letting $\mu_{j}^{*}=0$ for all $j \notin A\left(x^{*}\right)$ the desired relations (4.4) and (4.5) follow.
yields then a corresponding Lagrange multiplier result, where the Lagrange multipliers corresponding to equality constraints have unrestricted sign. In particular, if $x^{*}$ is a local minimum for the problem

$$
\begin{aligned}
\operatorname{minimize} & f(x) \\
\text { subject to } & e_{i}^{\prime} x=d_{i}, \quad i=1, \ldots, m \\
& a_{j}^{\prime} x \leq b_{j}, \quad j=1, \ldots, r
\end{aligned}
$$

---

then there exist $\lambda_{1}^{*}, \ldots, \lambda_{m}^{*}, \mu_{1}^{*}, \ldots, \mu_{r}^{*}$, such that

$$
\begin{gathered}
\nabla f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*} e_{i}+\sum_{j=1}^{r} \mu_{j}^{*} a_{j}=0 \\
\mu_{j}^{*} \geq 0 . \quad j=1, \ldots r \\
\mu_{j}^{*}=0 . \quad j \notin A\left(x^{*}\right)
\end{gathered}
$$

These conditions will hold even if $x^{*}$ is not regular. The Lagrange multipliers $\lambda_{i}^{*}$ are given by

$$
\lambda_{i}^{*}=w_{i}^{*}-v_{i}^{*}
$$

where $w_{i}^{*} \geq 0$ is a Lagrange multiplier corresponding to the constraint $e_{i}^{\prime} x-d_{i} \leq 0$ and $v_{i}^{*} \geq 0$ is a Lagrange multiplier corresponding to $-e_{i}^{\prime} x+$ $d_{i} \leq 0$.

# 3.4.1 Convex Cost Functions and Linear Constraints 

As one might expect based on the corresponding results of Chapter 2. if $f$ is convex, the necessary condition of Prop. 3.4.1 is also sufficient. Indeed, the condition $-\nabla f\left(x^{*}\right)=\sum_{j=1}^{r} \mu_{j}^{*} a_{j}$ implies that $\nabla f\left(x^{*}\right)^{\prime} y \geq 0$ for all the feasible directions at $x^{*}$, which are the vectors $y \neq 0$ with $a_{j}^{\prime} y \leq 0$ for all $j \in A\left(x^{*}\right)$. Thus if $f$ is convex, by Prop. 2.1.2(b) in Section 2.1. $x^{*}$ is a global minimum. In fact a sharper necessary and sufficient condition can be derived, involving minimization of a Lagrangian function. This minimization may involve any subset of the inequality constraints, while the remaining constraints are taken into account by the Lagrange multipliers [see Eq. (4.9) that follows]. This flexibility of assigning Lagrange multipliers to only some of the constraints, while dealing with the other constraints explicitly is often very useful.

Proposition 3.4.2: (Optimality Conditions for Convex Cost and Linear Constraints) Let $f: \Re^{n} \mapsto \Re$ convex and continuously differentiable, and let $J$ be a subset of the index set $\{1, \ldots, r\}$. Then $x^{*}$ is a global minimum for the problem

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } a_{j}^{\prime} x \leq b_{j}, \quad j=1, \ldots, r
\end{aligned}
$$

if and only if $x^{*}$ is feasible and there exist scalars $\mu_{j}^{*}, j \in J$, such that

$$
\mu_{j}^{*} \geq 0, \quad j \in J
$$

---

$$
\begin{gathered}
\mu_{j}^{*}=0, \quad \forall j \in J \text { with } j \notin A\left(x^{*}\right) \\
x^{*}=\arg \min _{\substack{a_{j}^{\prime} x \leq b_{j} \\
j \notin J}}\left\{f(x)+\sum_{j \in J} \mu_{j}^{*}\left(a_{j}^{\prime} x-b_{j}\right)\right\}
\end{gathered}
$$

Proof: Assume that $x^{*}$ is a global minimum. Then by Prop. 3.4.1, there exist nonnegative scalars $\mu_{1}^{*}, \ldots, \mu_{r}^{*}$ such that

$$
\mu_{j}^{*}\left(a_{j}^{\prime} x^{*}-b_{j}\right)=0, \quad j=1, \ldots, r
$$

and

$$
\nabla f\left(x^{*}\right)+\sum_{j=1}^{r} \mu_{j}^{*} a_{j}=0
$$

Using the convexity of $f$, the last relation implies that

$$
x^{*}=\arg \min _{x \in \Re^{n}}\left\{f(x)+\sum_{j=1}^{r} \mu_{j}^{*}\left(a_{j}^{\prime} x-b_{j}\right)\right\}
$$

Combining Eqs. (4.10) and (4.11), we obtain

$$
f\left(x^{*}\right)=\min _{x \in \Re^{n}}\left\{f(x)+\sum_{j=1}^{r} \mu_{j}^{*}\left(a_{j}^{\prime} x-b_{j}\right)\right\}
$$

Since $\mu_{j}^{*} \geq 0$, it follows that $\mu_{j}^{*}\left(a_{j}^{\prime} x-b_{j}\right) \leq 0$ if $a_{j}^{\prime} x-b_{j} \leq 0$, so Eq. (4.12) implies that

$$
\begin{aligned}
f\left(x^{*}\right) & \leq \min _{\substack{a_{j}^{\prime} x \leq b_{j} \\
j \notin J}}\left\{f(x)+\sum_{j=1}^{r} \mu_{j}^{*}\left(a_{j}^{\prime} x-b_{j}\right)\right\} \\
& \leq \min _{\substack{a_{j}^{\prime} x \leq b_{j} \\
j \notin J}}\left\{f(x)+\sum_{j \in J} \mu_{j}^{*}\left(a_{j}^{\prime} x-b_{j}\right)\right\}
\end{aligned}
$$

Equations (4.10) and (4.13) imply that $x^{*}$ attains the minimum in the right-hand side of the above relation as desired, thus proving the desired Eq. (4.9).

Conversely, assume that $x^{*}$ is feasible and there exist scalars $\mu_{j}^{*}, j \in$ $J$, satisfying conditions (4.7)-(4.9). Using Prop. 2.1.2(a) in Section 2.1, we

---

see that Eq. (4.9) implies

$$
\left(\nabla f\left(x^{*}\right)+\sum_{j \in J} \mu_{j}^{*} a_{j}\right)^{\prime}\left(x-x^{*}\right) \geq 0, \quad \forall x \text { with } a_{j}^{\prime} x \leq b_{j}, \forall j \notin J
$$

For all feasible vectors $x$ for the original problem $\min _{a_{j}^{\prime} x \leq b_{j}, j=1, \ldots, r} f(x)$, we have

$$
a_{j}^{\prime} x \leq b_{j}=a_{j}^{\prime} x^{*}, \quad \forall j \in A\left(x^{*}\right)
$$

Since $\mu_{j}^{*}=0$ if $j \in J$ and $j \notin A\left(x^{*}\right)$, we obtain

$$
\sum_{j \in J} \mu_{j}^{*} a_{j}^{\prime}\left(x-x^{*}\right) \leq 0
$$

and by using condition (4.14), we see that

$$
\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right) \geq 0
$$

for all feasible $x$. It follows using the convexity of $f$ and Prop. 2.1.2(b) in Section 2.1, that $x^{*}$ is a global minimum. Q.E.D.

Note from the preceding proof that given $x^{*}$, the vectors $\mu^{*}$ satisfying conditions (4.7)-(4.9) in the preceding proof are independent of the subset $J$; they are the vectors $\mu^{*}$ corresponding to $J=\{1, \ldots, r\}$, but they satisfy condition (4.9) for any subset $J$.

# Example 4.1 (Optimization Over a Simplex) 

Consider the case where $f$ is convex and the constraint'set is the simplex

$$
\left\{x \mid x \geq 0, \sum_{i=1}^{n} x_{i}=r\right\}
$$

where $r>0$ is a given scalar. We will use the optimality condition of Prop. 3.4.2 to recover the conditions obtained in Example 1.2 of Section 2.1.

By Prop. 3.4.2(a), if $x^{*}$ is a global minimum, there exists a scalar $\lambda^{*}$ such that

$$
x^{*}=\arg \min _{x \geq 0}\left\{f(x)+\lambda^{*}\left(r-\sum_{i=1}^{n} x_{i}\right)\right\}
$$

(We are applying here the extended form of the proposition that involves linear equality constraints as well; this is the same as the one involving inequality constraints exclusively, except that the sign of the Lagrange multipliers corresponding to the equality constraints is unrestricted.) By applying the

---

optimality conditions for an orthant constraint (cf. Example 1.1 in Section 2.1), we obtain from Eq. (4.15)

$$
\begin{gathered}
\frac{\partial f\left(x^{*}\right)}{\partial x_{i}} \geq \lambda^{*} ; \quad i=1, \ldots, n \\
\frac{\partial f\left(x^{*}\right)}{\partial x_{i}}=\lambda^{*} ; \quad \text { if } x_{i}^{*}>0
\end{gathered}
$$

We see therefore that all coordinates $x_{i}^{*}$ which are positive must have partial cost derivatives which are minimal and equal to the Lagrange multiplier $\lambda^{*}$.

By applying Prop. 3.4.2, we also see that the conditions (4.16) and (4.17) are sufficient for a feasible vector $x^{*}$ to be optimal.

# 3.4.2 Duality Theory: A Simple Form for Linear Constraints 

Given a problem with convex cost function and linear constraints, it is possible to construct another problem, called dual, that has the same optimal value and has as optimal solutions the Lagrange multipliers of the original. This is a very important relationship, which will be discussed extensively and generalized considerably in Chapter 5. Here, we develop some of the simpler aspects of duality theory, which follow almost immediately from the theory developed so far in this section.

Consider the problem

$$
\begin{aligned}
& \text { minimize } f(x) \\
& \text { subject to } \quad e_{i}^{\prime} x=d_{i}, \quad i=1, \ldots, m \\
& a_{j}^{\prime} x \leq b_{j}, \quad j=1, \ldots, r, \quad x \in X
\end{aligned}
$$

where $e_{i}, a_{j}$, and $d_{i}, b_{j}$ are given vectors and scalars, respectively, $f: \Re^{n} \rightarrow$ $\Re$ is a convex continuously differentiable function, and $X$ is a polyhedral set, that is, a set specified by a finite collection of linear equality and inequality constraints. We refer to problem (P) as the primal problem.

Define the Lagrangian function

$$
L(x, \lambda, \mu)=f(x)+\sum_{i=1}^{m} \lambda_{i}\left(e_{i}^{\prime} x-d_{i}\right)+\sum_{j=1}^{r} \mu_{j}\left(a_{j}^{\prime} x-b_{j}\right)
$$

Consider also the dual function $q: \Re^{m+r} \mapsto[-\infty, \infty)$ defined by

$$
q(\lambda, \mu)=\inf _{x \in X} L(x, \lambda, \mu)
$$

The dual problem is

$$
\begin{aligned}
& \text { maximize } q(\lambda, \mu) \\
& \text { subject to } \lambda \in \Re^{m}, \quad \mu \geq 0
\end{aligned}
$$

---

Note that if the polyhedral set $X$ is bounded, then it is also compact (every polyhedral set is closed since it is the intersection of closed subspaces), and the infimum in Eq. (4.19) is attained for all $(\lambda, \mu)$ by Weierstrass' theorem (Prop. A. 8 in Appendix A). Thus if $X$ is bounded, the dual function takes real values. In general, however, $q(\lambda, \mu)$ can take the value $-\infty$. Thus in effect, the constraint set of the dual problem (D) is the set

$$
Q=\{(\lambda, \mu) \mid \mu \geq 0, q(\lambda, \mu)>-\infty\}
$$

We have the following basic result:

# Proposition 3.4.3: (Duality Theorem) 

(a) If the primal problem ( P ) has an optimal solution, the dual problem (D) also has an optimal solution and the corresponding optimal values are equal.
(b) In order for $x^{*}$ to be an optimal primal solution and $\left(\lambda^{*}, \mu^{*}\right)$ to be an optimal dual solution, it is necessary and sufficient that $x^{*}$ is feasible for the primal problem, $\mu^{*} \geq 0$, and

$$
f\left(x^{*}\right)=L\left(x^{*}, \lambda^{*}, \mu^{*}\right)=\min _{x \in X} L\left(x, \lambda^{*}, \mu^{*}\right)
$$

Proof: (a) For all primal feasible $x$, and all $\lambda \in \Re^{m}$ and $\mu \geq 0$, we have $\lambda_{i}\left(e_{i}^{\prime} x-d_{i}\right)=0$ and $\mu_{j}^{\prime}\left(a_{j}^{\prime} x-b_{j}\right) \leq 0$, so that

$$
q(\lambda, \mu) \leq f(x)+\sum_{i=1}^{m} \lambda_{i}\left(e_{i}^{\prime} x-d_{i}\right)+\sum_{j=1}^{r} \mu_{j}\left(a_{j}^{\prime} x-b_{j}\right) \leq f(x)
$$

By taking the minimum of the right-hand side over all primal feasible $x$, we obtain

$$
q(\lambda, \mu) \leq f\left(x^{*}\right), \quad \forall \lambda \in \Re^{m}, \mu \geq 0
$$

where $x^{*}$ is a primal optimal solution. By Prop. 3.4.2, there exist $\lambda^{*} \in \Re^{m}$ and $\mu^{*} \geq 0$ such that $\mu_{j}^{*}\left(a_{j}^{\prime} x^{*}-b_{j}\right)=0$ for all $j$, and

$$
x^{*}=\arg \min _{x \in X} L\left(x, \lambda^{*}, \mu^{*}\right)
$$

so by using also the definition (4.19) of $q$, we have

$$
\begin{aligned}
q\left(\lambda^{*}, \mu^{*}\right) & =L\left(x^{*}, \lambda^{*}, \mu^{*}\right) \\
& =f\left(x^{*}\right)+\sum_{i=1}^{m} \lambda_{i}^{*}\left(e_{i}^{\prime} x^{*}-d_{i}\right)+\sum_{j=1}^{r} \mu_{j}^{*}\left(a_{j}^{\prime} x^{*}-b_{j}\right) \\
& =f\left(x^{*}\right)
\end{aligned}
$$

---

By combining Eqs. (4.22) and (4.23), we see that $\left(\lambda^{\star}, \mu^{\star}\right)$ is a dual optimal solution and that $q\left(\lambda^{\star}, \mu^{\star}\right)=f\left(x^{\star}\right)$.
(b) If $x^{\star}$ is primal optimal and $\left(\lambda^{\star}, \mu^{\star}\right)$ is dual optimal, by part (a) we obtain

$$
f\left(x^{\star}\right)=q\left(\lambda^{\star}, \mu^{\star}\right)
$$

which when combined with Eq. (4.21), yields

$$
f\left(x^{\star}\right)=L\left(x^{\star}, \lambda^{\star}, \mu^{\star}\right)=q\left(\lambda^{\star}, \mu^{\star}\right)=\min _{x \in X} L\left(x, \lambda^{\star}, \mu^{\star}\right)
$$

Conversely, the relation $f\left(x^{\star}\right)=\min _{x \in X} L\left(x, \lambda^{\star}, \mu^{\star}\right)$ can be written as $f\left(x^{\star}\right)=q\left(\lambda^{\star}, \mu^{\star}\right)$, and since $x^{\star}$ is primal feasible and $\mu^{\star} \geq 0$, Eq. (4.21) implies that $x^{\star}$ is primal optimal and $\left(\lambda^{\star}, \mu^{\star}\right)$ is dual optimal. Q.E.D.

The duality assertions of Prop. 3.4.3 were proved under hypotheses which are stronger than necessary. For example, the differentiability assumption on $f$ is not needed - it turns out that convexity is sufficient. Furthermore, the existence of a primal optimal solution is not needed for the conclusion of Prop. 3.4.3(a) - it turns out that finiteness of the infimum of the primal cost is sufficient. However, to prove these and other extensions of the duality theorem we need an approach which is not calculus-based like the one of this chapter but rather is based on the geometry of convex sets; see Chapter 5.

# Example 4.2 (The Dual of a Linear Program) 

Consider the linear program

$$
\begin{aligned}
& \operatorname{minimize} c^{\prime} x \\
& \text { subject to } e_{i}^{\prime} x=d_{i}, \quad i=1, \ldots, m, \quad x \geq 0
\end{aligned}
$$

where $c$ and $e_{i}$ are given vectors in $\Re^{n}$, and $d_{i}$ are given scalars. We consider the dual function

$$
q(\lambda)=\inf _{x \geq 0}\left\{\sum_{j=1}^{n}\left(c_{j}-\sum_{i=1}^{m} \lambda_{i} e_{i j}\right) x_{j}+\sum_{i=1}^{m} \lambda_{i} d_{i}\right\}
$$

where $e_{i j}$ is the $j$ th coordinate of the vector $e_{i}$. It is seen that if $c_{j}-$ $\sum_{i=1}^{m} \lambda_{i} e_{i j} \geq 0$ for all $j$, the infimum above is attained for $x=0$, and we have $q(\lambda)=\sum_{i=1}^{m} \lambda_{i} d_{i}$. On the other hand, if $c_{j}-\sum_{i=1}^{m} \lambda_{i} e_{i j}<0$ for some $j$, we can make the expression in braces above arbitrarily small by taking $x_{j}$ sufficiently large, so that we have $q(\lambda)=-\infty$ in this case. Thus, the dual problem is

$$
\begin{aligned}
& \operatorname{maximize} \sum_{i=1}^{m} \lambda_{i} d_{i} \\
& \text { subject to } \sum_{i=1}^{m} \lambda_{i} e_{i j} \leq c_{j}, \quad j=1, \ldots, n
\end{aligned}
$$

---

By the optimality conditions of Prop. 3.4.3, $\left(x^{*}, \lambda^{*}\right)$ is a primal and dual optimal solution pair if and only if $x^{*}$ is primal feasible, $\lambda^{*}$ is dual feasible, and the Lagrangian optimality condition

$$
x^{*}=\arg \min _{x \geq 0}\left\{\left(c-\sum_{i=1}^{m} \lambda_{i}^{*} e_{i}\right)^{\prime} x+\sum_{i=1}^{m} \lambda_{i}^{*} d_{i}\right\}
$$

holds [cf. Eq. (4.20)]. Given the dual feasibility condition on $\lambda^{*}$, the Lagrangian optimality condition is equivalent to the following two relations

$$
\begin{aligned}
& x_{j}^{*}>0 \quad \Rightarrow \quad \sum_{i=1}^{m} \lambda_{i}^{*} e_{i j}=c_{j}, \quad j=1, \ldots, n, \\
& \sum_{i=1}^{m} \lambda_{i}^{*} e_{i j}<c_{j} \quad \Rightarrow \quad x_{j}^{*}=0, \quad j=1, \ldots, n,
\end{aligned}
$$

known as the complementary slackness conditions.
Let us consider now the dual of the dual problem (DLP). We first convert this program into the equivalent minimization problem

$$
\begin{aligned}
& \operatorname{minimize} \sum_{i=1}^{m}\left(-d_{i}\right) \lambda_{i} \\
& \text { subject to } \sum_{i=1}^{m} \lambda_{i} e_{i j} \leq c_{j}, \quad j=1, \ldots, n
\end{aligned}
$$

Assigning a Lagrange multiplier $x_{j}$ to the $j$ th inequality constraint, the dual function of this problem is given by

$$
\begin{aligned}
p(x) & =\min _{\lambda \in \Re^{m}}\left\{\sum_{i=1}^{m}\left(\sum_{j=1}^{n} e_{i j} x_{j}-d_{i}\right) \lambda_{i}-\sum_{j=1}^{n} c_{j} x_{j}\right\} \\
& = \begin{cases}-c^{\prime} x, & \text { if } e_{i}^{\prime} x=d_{i}, \quad i=1, \ldots, m \\
-\infty, & \text { otherwise. }\end{cases}
\end{aligned}
$$

The corresponding dual problem is

$$
\begin{aligned}
& \text { maximize } p(x) \\
& \text { subject to } x \geq 0
\end{aligned}
$$

or equivalently

$$
\begin{aligned}
& \text { minimize } c^{\prime} x \\
& \text { subject to } e_{i}^{\prime} x=d_{i}, \quad i=1, \ldots, m, \quad x \geq 0
\end{aligned}
$$

which is identical to the primal problem (LP). We have thus shown that the duality is symmetric, that is, the dual of the dual linear program (DLP) is the primal problem $(L P)$.

---

The pair of primal and dual linear programs (LP) and (DLP) can also be written compactly in terms of the $m \times n$ matrix $E$ having rows $e_{1}^{\prime}, \ldots, e_{m}^{\prime}$, and the vector $d$ having coordinates $d_{1}, \ldots, d_{m}$;

$$
\min _{E x=d, x \geq 0} c^{\prime} x \quad \Longleftrightarrow \quad \max _{E^{\prime} \lambda \leq c} d^{\prime} \lambda
$$

Other linear programming duality relations that can be verified by the reader include the following:

$$
\begin{aligned}
\min _{A^{\prime} x \geq b} c^{\prime} x & \Longleftrightarrow \max _{A \mu=c, \mu \geq 0} b^{\prime} \mu \\
\min _{A^{\prime} x \geq b, x \geq 0} c^{\prime} x & \Longleftrightarrow \max _{A \mu \leq c, \mu \geq 0} b^{\prime} \mu
\end{aligned}
$$

# Example 4.3 (The Dual of a Quadratic Program) 

Consider the quadratic programming problem

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2} x^{\prime} Q x+c^{\prime} x \\
& \text { subject to } A x \leq b
\end{aligned}
$$

where $Q$ is a given $n \times n$ positive definite symmetric matrix, $A$ is a given $r \times n$ matrix, and $b \in \Re^{r}$ and $c \in \Re^{n}$ are given vectors. The dual function is

$$
q(\mu)=\inf _{x \in \Re^{n}}\left\{\frac{1}{2} x^{\prime} Q x+c^{\prime} x+\mu^{\prime}(A x-b)\right\}
$$

The infimum is attained for $x=-Q^{-1}\left(c+A^{\prime} \mu\right)$, and, after substitution of this expression in the preceding relation for $q$, a straightforward calculation yields

$$
q(\mu)=-\frac{1}{2} \mu^{\prime} A Q^{-1} A^{\prime} \mu-\mu^{\prime}\left(b+A Q^{-1} c\right)-\frac{1}{2} c^{\prime} Q^{-1} c
$$

The dual problem, after a sign change that converts it to a minimization problem, can be written as

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2} \mu^{\prime} P \mu+t^{\prime} \mu \\
& \text { subject to } \mu \geq 0
\end{aligned}
$$

where

$$
P=A Q^{-1} A^{\prime}, \quad t=b+A Q^{-1} c
$$

If $\mu^{*}$ is any dual optimal solution, then the optimal solution of the primal problem (QP) is

$$
x^{*}=-Q^{-1}\left(c+A^{\prime} \mu^{*}\right)
$$

Note that the dual problem is also a quadratic programming problem, but it has simpler constraints than the primal. Furthermore, if the row dimension $r$ of $A$ is smaller than its column dimension $n$, the dual problem is defined on a space of smaller dimension than the primal.

---

# 3.4.3 Relation Between Primal and Dual Functions Sensitivity 

We can use duality theory to obtain sharper sensitivity results under the convexity assumptions of this section. To this end, let us define the primal function $p: \Re^{m+n} \mapsto(-\infty, \infty]$ for problem (P) by

$$
p(u, v)=\inf _{\substack{x \in X \\ c_{x}^{\prime} x=d_{i}+u_{i-1}=1, \ldots, m \\ a_{j}^{\prime} x \leq b_{j}+v_{j}, j=1, \ldots, r}} f(x)
$$

where $u=\left(u_{1}, \ldots, u_{m}\right)$ and $v=\left(v_{1}, \ldots, v_{r}\right)$ are the vectors that specify the perturbations of the right-hand sides of the constraints. Note that there may exist some ( $u, v$ ) such that the minimization problem above has no feasible solution. For such $(u, v)$, we define $p(u, v)=\infty$.

The preceding definition of the primal function $p$ differs from the ones of Sections 3.2 and 3.3 in that $p$ is defined without reference to a particular local minimum. Under the convexity assumptions of this subsection, the primal function $p$, the dual function $q$, and the corresponding sets

$$
P=\{(u, v) \mid p(u, v)<\infty\}
$$

and

$$
Q=\{(\lambda, \mu) \mid \mu \geq 0, q(\lambda, \mu)>-\infty\}
$$

have certain convexity properties as shown in the following proposition.

Proposition 3.4.4: The sets $P$ and $Q$ of Eqs. (4.28) and (4.29) are convex. Furthermore, the primal function $p$ is convex over $P$ and the dual function $q$ is concave over $Q$.

Proof: To simplify notation, we assume in this proof that there are no equality constraints, so that $p$ is a function of $v$ only and $q$ is a function of $\mu$ only. To see that $P$ is convex, let $v^{1} \in P$ and $v^{2} \in P$ and let $x^{1} \in X$ and $x^{2} \in X$ be such that $a_{j}^{\prime} x^{1} \leq b_{j}+v_{j}^{1}$ and $a_{j}^{\prime} x^{2} \leq b_{j}+v_{j}^{2}$ for all $j$. Then we have for all $j$ and $\alpha \in[0,1]$,

$$
a_{j}^{\prime}\left(\alpha x^{1}+(1-\alpha) x^{2}\right) \leq b_{j}+\alpha v_{j}^{1}+(1-\alpha) v_{j}^{2}
$$

so it follows that $\alpha v^{1}+(1-\alpha) v^{2} \in P$ and $P$ is convex. Furthermore, if $\left\{x_{1}^{k}\right\} \subset X$ and $\left\{x_{2}^{k}\right\} \subset X$ are sequences such that

$$
\begin{gathered}
a_{j}^{\prime} x_{1}^{k} \leq b_{j}+v_{j}^{1}, \quad a_{j}^{\prime} x_{2}^{k} \leq b_{j}+v_{j}^{2}, \quad \forall j=1, \ldots, r, k=0,1, \ldots \\
f\left(x_{1}^{k}\right) \rightarrow p\left(v^{1}\right), \quad f\left(x_{2}^{k}\right) \rightarrow p\left(v^{2}\right)
\end{gathered}
$$

---

then using the convexity of $f$, we have

$$
\begin{aligned}
p\left(\alpha v^{1}+(1-\alpha) v^{2}\right) & \leq \liminf _{k \rightarrow \infty} f\left(\alpha x_{1}^{k}+(1-\alpha) x_{2}^{k}\right) \\
& \leq \lim _{k \rightarrow \infty}\left\{\alpha f\left(x_{1}^{k}\right)+(1-\alpha) f\left(x_{2}^{k}\right)\right\} \\
& =\alpha p\left(v^{1}\right)+(1-\alpha) p\left(v^{2}\right)
\end{aligned}
$$

proving the convexity of $p$.
For any $x, \mu, \bar{\mu}$, and $\alpha \in[0,1]$, we have

$$
L(x, \alpha \mu+(1-\alpha) \bar{\mu})=\alpha L(x, \mu)+(1-\alpha) L(x, \bar{\mu})
$$

Taking infimum over $x \in X$, we obtain

$$
\inf _{x \in X} L(x, \alpha \mu+(1-\alpha) \bar{\mu}) \geq \alpha \inf _{x \in X} L(x, \mu)+(1-\alpha) \inf _{x \in X} L(x, \bar{\mu})
$$

or

$$
q(\alpha \mu+(1-\alpha) \bar{\mu}) \geq \alpha q(\mu)+(1-\alpha) q(\bar{\mu})
$$

This implies that if $\mu$ and $\bar{\mu}$ belong to $Q$, the same is true for $\alpha \mu+(1-\alpha) \bar{\mu}$ and hence $Q$ is convex. Equation (4.30) also implies that $q$ is concave over Q. Q.E.D.

We next show how the primal function is related to the dual function.

Proposition 3.4.5: Assume that for all $(u, v) \in P$, the problem

$$
\begin{aligned}
& \text { minimize } f(x) \\
& \text { subject to } \quad e_{i}^{\prime} x=d_{i}+v_{i}, \quad i=1, \ldots, m \\
& a_{j}^{\prime} x \leq b_{j}+v_{j}, \quad j=1, \ldots, r, \quad r \in X
\end{aligned}
$$

has an optimal solution. Then the following hold for the primal function

$$
p(u, v)=\inf _{\substack{x \in X \\ e_{i}^{\prime} x=d_{i}+u_{i}, i=1, \ldots, m \\ a_{j}^{\prime} x \leq b_{j}+v_{j}, j=1, \ldots, r}}
$$

and the dual function

$$
q(\lambda, \mu)=\inf _{x \in X}\left\{f(x)+\sum_{i=1}^{m} \lambda_{i}\left(e_{i}^{\prime} x-d_{i}\right)+\sum_{j=1}^{r} \mu_{j}\left(a_{j}^{\prime} x-b_{j}\right)\right\}
$$

---

(a) For all $(u, v) \in P$ we have

$$
p(u, v)=\max _{(\lambda, \mu) \in Q}\left\{q(\lambda, \mu)-\lambda^{\prime} u-\mu^{\prime} v\right\}
$$

(b) For all $(\lambda, \mu) \in Q$ we have

$$
q(\lambda, \mu)=\inf _{(u, v) \in P}\left\{p(u, v)+\lambda^{\prime} u+\mu^{\prime} v\right\}
$$

Furthermore, ( $u, v$ ) attains the minimum in Eq. (4.33) if and only if $(\lambda, \mu)$ attains the maximum in Eq. (4.32).
(c) For a given $(u, v) \in P,(\lambda, \mu)$ is an optimal dual solution for problem (4.31) if and only if $-(\lambda, \mu)$ is a subgradient of the primal function $p$ at $(u, v)$.
(d) For a given $(u, v) \in P,(\lambda, \mu)$ is an optimal dual solution for problem (4.31) if and only if $(u, v)$ is a subgradient of the dual function $q$ at $(\lambda, \mu)$.

Proof: To simplify notation, we assume that there are no equality constraints, so that $p$ is a function of $r$ only and $q$ is a function of $\mu$ only.
(a) For all $v \in P$, the dual function of problem (4.31) is

$$
q_{v}(\mu)=\inf _{x \in X}\left\{f(x)+\sum_{j=1}^{r} \mu_{j}\left(a_{j}^{\prime} x-b_{j}-v_{j}\right)\right\}=q(\mu)-\mu^{\prime} v
$$

Since by assumption, problem (4.31) has an optimal solution, Prop. 3.4.3(a) shows that $p(v)$ is equal to the maximal value of $q_{v}(\mu)$ over $\mu \in Q$. This, together with the above equation, proves the result.
(b) For every $\mu \geq 0$, we have

$$
\begin{aligned}
q(\mu) & =\inf _{x \in X}\left\{f(x)+\sum_{j=1}^{r} \mu_{j}\left(a_{j}^{\prime} x-b_{j}\right)\right\} \\
& =\inf _{v \in P, x \in X, a_{j}^{\prime} x-b_{j} \leq v_{j}, j=1, \ldots, r}\left\{f(x)+\sum_{j=1}^{r} \mu_{j}\left(a_{j}^{\prime} x-b_{j}\right)\right\}
\end{aligned}
$$

Since $\mu \geq 0$ and $P$ contains all the vectors of the form $v=\left(a_{1}^{\prime} x-\right.$

---

$b_{1}, \ldots, a_{r}^{\prime} x-b_{r}$ ) with $x \in X$, this relation yields

$$
\begin{aligned}
q(\mu) & =\inf _{v \in P, x \in X, a_{j}^{\prime} x-b_{j} \leq v_{j}, j=1, \ldots, r}\left\{f(x)+\sum_{j=1}^{r} \mu_{j} v_{j}\right\} \\
& =\inf _{v \in P} \inf _{x \in X, a_{j}^{\prime} x-b_{j} \leq v_{j}, j=1, \ldots, r}\left\{f(x)+\sum_{j=1}^{r} \mu_{j} v_{j}\right\}
\end{aligned}
$$

and finally

$$
q(\mu)=\inf _{v \in P}\left\{p(v)+\mu^{\prime} v\right\}, \quad \forall \mu \geq 0
$$

Furthermore, we have that $v$ attains the minimum in the above equation if and only if $v \in P$ and $q(\mu)=p(v)+\mu^{\prime} v$, which in view of Eq. (4.32), is true if and only if $\mu$ attains the maximum of $q(\bar{\mu})-\bar{\mu}^{\prime} v$ over $\bar{\mu} \in Q$.
(c) By part (a), for a given $v \in P, \mu$ is an optimal dual solution of problem (4.31) if and only if

$$
p(v)=q(\mu)-\mu^{\prime} v
$$

which together with Eq. (4.33) yields

$$
p(v) \leq p(\bar{v})+\mu^{\prime} \bar{v}-\mu^{\prime} v, \quad \forall \bar{v} \in P
$$

This is equivalent to $-\mu$ being a subgradient of $p$ at $v$ [cf. Eq. (B.28) in Appendix B].
(d) By parts (a) and (b), for a given $v \in P, \mu$ is an optimal dual solution of problem (4.31) if and only if

$$
q(\mu)=p(v)+\mu^{\prime} v
$$

which together with Eq. (4.32) yields

$$
q(\mu) \geq q(\bar{\mu})-\bar{\mu}^{\prime} v+\mu^{\prime} v, \quad \forall \bar{\mu} \in Q
$$

This is equivalent to $v$ being a subgradient of $q$ at $\mu$. Q.E.D.
The relation between the primal and dual functions given in Eqs. (4.32) and (4.33) is called a conjugacy relation and will be studied in greater detail in Section 5.4. Part (c) of Prop. 3.4.5 generalizes the sensitivity results of Sections 3.2 and 3.3. It asserts in particular that if $p$ is differentiable at $(0,0)$, there exists a unique optimal dual solution $(\lambda, \mu)$; this solution can be viewed as a Lagrange multiplier and satisfies

$$
\nabla_{u} p(0,0)=-\lambda, \quad \nabla_{v} p(0,0)=-\mu
$$

The relation between the primal and dual functions is sometimes useful when we want to evaluate $p(u, v)$ for a range of values of $(u, v)$. We can instead evaluate $q(\lambda, \mu)$ for an appropriate range of values of $(\lambda, \mu)$ and obtain the corresponding vectors $(u, v)$ as subgradients of $q$ in accordance with Prop. 3.4.5(d). The following example illustrates this situation.

---

# Example 4.4 (Linear Programming Sensitivity) 

Consider a linear program of the form

$$
\begin{aligned}
& \operatorname{minimize} c^{\prime} x \\
& \text { subject to } x \in X, \quad a^{\prime} x \leq b
\end{aligned}
$$

where $X$ is a bounded polyhedral set in $\Re^{\prime \prime}$, and $a^{\prime} x \leq b$ is a single inequality constraint. Suppose we are interested in sensitivity analysis as the right-hand side of this constraint is perturbed. The corresponding primal function is

$$
p(v)=\min _{\substack{x \in X \\ a^{\prime} x \leq b+v}} c^{\prime} x
$$

Because the polyhedron $X$ is bounded, it is also compact, and the minimization problem above has an optimal solution if and only if it has a feasible solution, that is, $v \in P$ (Weierstrass' theorem, given as Prop. A. 8 in Appendix A). Furthermore, Prop. 3.4.5 applies.

The dual function

$$
q(\mu)=\min _{x \in X}\left\{c^{\prime} x+\mu\left(a^{\prime} x-b\right)\right\}
$$

is a scalar piecewise linear function, since the minimum in the above equation is attained at one of the finitely many extreme points of $X$ (Prop. B. 21 in Appendix B). In view of the relation between $p$ and $q$ expressed by Prop. 3.4.5, the primal function $p$ is also piecewise linear. Furthermore, it can be shown that the break points of $p$ correspond to the slopes of the linear segments of $q$, and the break points of $q$ correspond to the slopes of the linear segments of $q$.

Suppose now that for purposes of sensitivity analysis we want to solve the problem

$$
\begin{aligned}
& \operatorname{minimize} c^{\prime} x \\
& \text { subject to } x \in X, \quad a^{\prime} x \leq b+v
\end{aligned}
$$

for a range $\left[\beta_{1}, \beta_{2}\right]$ of values of $v$. This can also be done by solving instead the possibly easier problem

$$
\begin{aligned}
& \operatorname{minimize} c^{\prime} x+\mu\left(a^{\prime} x-b\right) \\
& \text { subject to } x \in X
\end{aligned}
$$

for a range of values $\left[\gamma_{1}, \gamma_{2}\right]$ of $\mu$. The value $\gamma_{1}$ must be such that the right slope of $q$ at $\gamma_{1}$ is at least $\beta_{2}$, and the value $\gamma_{2}$ must be such that the left slope of $q$ at $\gamma_{2}$ is at most $\beta_{1}$.

Let us also denote by $X_{v}$ and $X_{\mu}$ the sets of optimal solutions of problems (4.35) and (4.36), respectively. Then it can be seen that if $\mu$ is a subgradient of $p$ at $v$, we always have

$$
X_{\mu} \subset X_{v}
$$

---

and we have $X_{\mu}=X_{v}$ if and only if, in addition, $\mu$ is not a slope of a linear segment of the primal function $p$, or equivalently, $\mu$ is not one of the break points of the dual function $q$.

# E XERCISES 

## 4.1

A company has available for sale a quantity of $Q$ units of a certain product that will be sold in $n$ market outlets. For each $i=1, \ldots, n$, the quantity $d_{i}$, which is demanded at outlet $i$, and the price of sale $p_{i}$ are known. The company wishes to determine the quantities $s_{i}^{*}$ with $0 \leq s_{i}^{*} \leq d_{i}$ to be sold at each outlet that maximize the revenue $\sum_{i=1}^{n} p_{i} s_{i}$ from the sale. Assuming that $d_{i}>0, p_{i}>0$, and $\sum_{i=1}^{n} d_{i} \geq Q$, show that there exists a cutoff price level $y$ such that for each $i$ if $p_{i}>y$, then $s_{i}^{*}=d_{i}$ and if $p_{i}<y$, then $s_{i}^{*}=0$. What happens if $p_{i}=y$ ? Describe a procedure for obtaining $s_{i}^{*}$.

## 4.2

Consider a problem of finding the best way to place bets totalling $A$ dollars $(A>0)$ in a race involving $n$ horses. Assume that we know the probability $p_{i}$, that the $i$ th horse wins, and the amount $s_{i}>0$ that the rest of the public is betting on the $i$ th horse. The track keeps a proportion $1-c$ of the total amount $(0<1-c<1)$ and distributes the rest among the public in proportion to the amounts bet on the winning horse. Thus, if we bet $x_{i}$ on the $i$ th horse, we receive

$$
c\left(A+\sum_{i=1}^{n} s_{i}\right) \frac{x_{i}}{s_{i}+x_{i}}
$$

if the $i$ th horse wins. The problem is to find $x_{i}^{*}, i=1, \ldots, n$, which maximize the expected net return

$$
c\left(A+\sum_{i=1}^{n} s_{i}\right)\left\{\sum_{i=1}^{n} \frac{p_{i} x_{i}}{s_{i}+x_{i}}\right\}-A
$$

or, equivalently,

$$
\sum_{i=1}^{n} \frac{p_{i} x_{i}}{s_{i}+x_{i}}
$$

subject to $\sum_{i=1}^{n} x_{i}=A, x_{i} \geq 0, i=1, \ldots, n$.

---

Assume that

$$
\frac{p_{1}}{s_{1}}>\frac{p_{2}}{s_{2}}>\cdots>\frac{p_{n}}{s_{n}}
$$

and show that there exists a scalar $\lambda^{*}$ such that the optimal solution is

$$
x_{i}^{*}= \begin{cases}\sqrt{\frac{s_{i} p_{i}}{\lambda^{*}}}-s_{i} & \text { for } i=1, \ldots, m^{*} \\ 0 & \text { for } i=m^{*}+1, \ldots, n\end{cases}
$$

where $m^{*}$ is the largest index $m$ for which $\frac{p_{m}}{s_{m}} \geq \lambda^{*}$.

# 4.3 

Verify the linear programming duality relations

$$
\begin{aligned}
& \min _{A^{\prime} x \geq b} c^{\prime} x \quad \Longleftrightarrow \quad \max _{A \mu=c, \mu \geq 0} b^{\prime} \mu \\
& \min _{A^{\prime} x \geq b, x \geq 0} c^{\prime} x \quad \Longleftrightarrow \quad \max _{A \mu \leq c, \mu \geq 0} b^{\prime} \mu
\end{aligned}
$$

show that they are symmetric, and derive the corresponding complementary slackness conditions [cf. Eqs. (4.25) and (4.26)].

### 4.4 (Duality for Transportation Problems)

Suppose that a quantity of a certain material must be shipped from $m$ supply points to $n$ demand points so as to minimize the total transportation cost. The supply at point $i$ is denoted $\alpha_{i}$ and the demand at point $j$ is denoted $\beta_{j}$. The unit transportation cost from $i$ to $j$ is $a_{i j}$. Letting $x_{i j}$ denote the quantity shipped from supply point $i$ to demand point $j$, the problem is

$$
\begin{aligned}
& \operatorname{minimize} \sum_{i, j} a_{i j} x_{i j} \\
& \text { subject to } \\
& \sum_{j=1}^{n} x_{i j}=\alpha_{i}, \quad \forall i=1, \ldots, m \\
& \sum_{i=1}^{m} x_{i j}=\beta_{j}, \quad \forall j=1, \ldots, n \\
& 0 \leq x_{i j}, \quad \forall i, j
\end{aligned}
$$

where $\alpha_{i}$ and $\beta_{j}$ are positive scalars, which for feasibility must satisfy

$$
\sum_{i=1}^{m} \alpha_{i}=\sum_{j=1}^{n} \beta_{j}
$$

---

(a) Assign Lagrange multipliers to the equality constraints and derive the corresponding dual problem.
(b) Introduce a price $p_{j}$ that demand point $j$ will pay per unit delivered at $j$. Show that if $x^{*}$ is an optimal solution of the transportation problem, there is a set of prices $\left\{p_{j}^{*} \mid j=1, \ldots, n\right\}$ such that if $x_{i j}^{*}>0$ then $j$ offers maximum net profit for $i$, that is,

$$
p_{j}^{*}-a_{i j}=\max _{k=1, \ldots, n}\left\{p_{k}^{*}-a_{i k}\right\}
$$

(c) Relate prices to dual feasible solutions, and show that every dual optimal solution has a property of the type described in (b) above.

# 3.5 NOTES AND SOURCES 

Section 3.1: The material on Lagrange multipliers is classical. For further elaboration of the treatment of Lagrange multipliers using augmented Lagrangians, see [Hes75]. For an alternative proof of the Lagrange multiplier theorem, see [Lue84]. For a survey of Lagrange multipliers and optimality conditions for more general problems than the ones considered here, including nonconvex nondifferentiable problems, see [Roc93].
Section 3.2: For a textbook treatment of sensitivity analysis, see [Fia78]. For sensitivity analysis under assumptions that are weaker than the second order sufficiency conditions that we have assumed here, see [DoT86], [Rob87], [GaJ88], [Sha88], [AuC90], [Kyp90], [Bon92], [KiR92], and [Iof94].
Section 3.3: Lagrange multiplier results for inequality constraints became available considerably later than their equality constraint counterparts. Important early works are those of Karush [Kar39] (an unpublished MS thesis), John [Joh48], and Kuhn and Tucker [KuT51]. The survey [Kuh76] gives a historical account of the development of the subject. There has been considerable effort to derive optimality conditions under assumptions that are weaker than regularity. Such conditions are generically called constraint qualifications (see Exercise 3.5). Important examples are those of [DuM65], [MaF67], and [Gui69]. For textbook treatments see [Man69] and [Hes75]. There has been much subsequent work on the subject, some of which addresses nondifferentiable problems, [GoT71], [BGS72], [GoT72], [Ben80], [BeZ82], [Cla83], [DeV85], [Mor88], and [Roc93].
Section 3.4: Duality theory has its origins in the work of von Neuman on zero sum games. The proof of linear programming duality was given by Gale, Kuhn, and Tucker [GKT51]. Most of the research on duality has been done under the weaker assumptions and the geometrical framework of Chapter 5; see the references in Section 5.5.

---
