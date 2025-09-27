# Optimization over a Convex Set 

## Contents

2.1. Optimality Conditions ..... p. 174
2.2. Feasible Directions and the Conditional Gradient Method ..... p. 191
2.3. Gradient Projection Methods ..... p. 203
2.4. Two-Metric Projection Methods ..... p. 224
2.5. Manifold Suboptimization - Quadratic Programming ..... p. 230
2.6. Affine Scaling for Linear Programming ..... p. 238
2.7. Block Coordinate Descent Methods* ..... p. 246
2.8. Notes and Sources ..... p. 250

---

In this chapter we consider the constrained optimization problem

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } x \in X
\end{aligned}
$$

where we assume throughout that:
(a) $X$ is a nonempty, convex, and closed subset of $\Re^{n}$.
(b) $f$ is continuously differentiable over $X$.

The problem of this chapter generalizes the unconstrained optimization problem of the preceding chapter, where $X=\Re^{n}$. We will see that the main algorithmic ideas for solving the unconstrained and the constrained problems are quite similar.

Usually the set $X$ has structure specified by equations and inequalities. If we take into account this structure, some new ideas, based on Lagrange multipliers and duality theory, come into play. These ideas will not be discussed in the present chapter, but they will be the subject of subsequent chapters.

Similar to the preceding chapter, the methods of this chapter are based on iterative descent along suitable directions. However, these directions must have the additional property that they maintain feasibility of the iterates. Such directions are called feasible, and as will be seen in Sections 2.2-2.7, they are usually obtained by solving certain optimization subproblems. As a general rule, these direction finding subproblems tend to be simpler when the constraint set $X$ is a polyhedron, and for this reason the methods of this chapter tend to be most suitable for linear constraints.

# 2.1 OPTIMALITY CONDITIONS 

In this section we consider the main necessary and sufficient optimality conditions for our problem, and we provide some examples of their application.

## Local and Global Minima

Throughout this book, a vector satisfying the constraints of a given problem will be called a feasible vector for that problem. A vector $x^{*} \in X$ is a local minimum of $f$ over the set $X$ if it is no worse than its feasible neighbors; that is, if there exists an $\epsilon>0$ such that

$$
f\left(x^{*}\right) \leq f(x), \quad \forall x \in X \text { with }\left\|x-x^{*}\right\|<\epsilon
$$

---

A vector $x^{*} \in X$ is a global minimum of $f$ over the set $X$ if it is no worse than all other feasible vectors, that is,

$$
f\left(x^{*}\right) \leq f(x), \quad \forall x \in X
$$

The local or global minimum $x^{*}$ is said to be strict if the corresponding inequality above is strict for $x \neq x^{*}$. Figure 2.1.1 illustrates these definitions. Local and global maxima of $f$ over $X$ are similarly defined; they are the local and global minima of the function $-f$ over $X$.

If both $f$ and $X$ are convex, a local minimum is also global. In particular, we have the following proposition, which is proved as Prop. B. 10 in Appendix B.

Proposition 2.1.1: If $f$ is a convex function, then a local minimum of $f$ over $X$ is a global minimum. If in addition $f$ is strictly convex over $X$, then there exists at most one global minimum of $f$ over $X$.
![[ch_2_p3_img1.jpeg]]

Figure 2.1.1. Local and global minima of $f$ over $X$.

# Necessary and Sufficient Conditions for Optimality 

As in unconstrained optimization, we expect that at a local minimum $x^{*}$, the first order variation $\nabla f\left(x^{*}\right)^{\prime} \Delta x$ due to a small feasible variation $\Delta x$ is nonnegative. Because $X$ is convex, the feasible variations are of the form

---

$\Delta x=x-x^{*}$, where $x \in X$. Thus the condition $\nabla f\left(x^{*}\right)^{\prime} \Delta x \geq 0$ translates into the necessary condition $\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right) \geq 0$ for all $x \in X$. The following proposition formally proves this condition and also shows that it is a sufficient condition for optimality when $f$ is convex.

# Proposition 2.1.2: (Optimality Condition) 

(a) If $x^{*}$ is a local minimum of $f$ over $X$, then

$$
\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right) \geq 0, \quad \forall x \in X
$$

(b) If $f$ is convex over $X$, then the condition of part (a) is also sufficient for $x^{*}$ to minimize $f$ over $X$.

Proof: (a) Suppose that $\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right)<0$ for some $x \in X$. By the [[Mean Value Theorem]] (Prop. A. 22 of Appendix A), for every $\epsilon>0$ there exists an $s \in[0,1]$ such that

$$
f\left(x^{*}+\epsilon\left(x-x^{*}\right)\right)=f\left(x^{*}\right)+\epsilon \nabla f\left(x^{*}+s \epsilon\left(x-x^{*}\right)\right)^{\prime}\left(x-x^{*}\right)
$$

Since $\nabla f$ is continuous, we have for all sufficiently small $\epsilon>0, \nabla f\left(x^{*}+\right.$ $\left.s \epsilon\left(x-x^{*}\right)\right)^{\prime}\left(x-x^{*}\right)<0$ and it follows that $f\left(x^{*}+\epsilon\left(x-x^{*}\right)\right)<f\left(x^{*}\right)$. The vector $x^{*}+\epsilon\left(x-x^{*}\right)$ is feasible for all $\epsilon \in[0,1]$ because $X$ is convex, so we obtain a contradiction of the local optimality of $x^{*}$.
(b) Using the convexity of $f$ and [[Mathematics/Optimization/Nonlinear Optimization/Textbook/Appendix#Characterizations of Differentiable Convex Functions|Prop. B. 3]] of Appendix B, we have $f(x) \geq$ $f\left(x^{*}\right)+\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right)$ for every $x \in X$. If the condition $\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right) \geq$ 0 holds for all $x \in X$, we obtain $f(x) \geq f\left(x^{*}\right)$ and therefore, $x^{*}$ minimizes $f$ over $X$. Q.E.D.
![[Mathematics/Optimization/Nonlinear Optimization/Textbook/Images/ch_2_p4_img2.jpeg]]

Figure 2.1.2. Geometric interpretation of the optimality condition of Prop. 2.1.2. At a local minimum $x^{*}$, the gradient $\nabla f\left(x^{*}\right)$ makes an angle less than or equal to 90 degrees with all feasible variations $x-x^{*}, x \in X$.

---

Figure 2.1.2 interprets the optimality condition (1.1) geometrically. Figure 2.1.3 shows how the condition can fail in the absence of convexity of $X$.

A vector $x^{*}$ satisfying the optimality condition (1.1) is referred to as a stationary point. In the absence of convexity of $f$, this condition may also be satisfied by local maxima and other points. To see this, note that if $X=\Re^{n}$ or if $x^{*}$ is an interior point of $X$, the condition (1.1) reduces to the stationarity condition $\nabla f\left(x^{*}\right)=0$ of unconstrained optimization. We now illustrate the condition (1.1) in some examples.
![[Mathematics/Optimization/Nonlinear Optimization/Textbook/Images/ch_2_p5_img3.jpeg]]

Figure 2.1.3. Illustration of how the necessary optimality condition may fail when $X$ is not convex. Here $x^{*}$ is a local minimum but we have $\nabla f\left(x^{*}\right)^{\prime}(x-$ $\left.x^{*}\right)<0$ for the feasible vector $x$ shown.

# Example 1.1 (Optimization Subject to Bounds on the Variables) 

Consider a positive orthant constraint

$$
X=\{x \mid x \geq 0\}
$$

Then the necessary condition (1.1) for $x^{*}=\left(x_{1}^{*}, \ldots, x_{n}^{*}\right)$ to be a local minimum is

$$
\sum_{i=1}^{n} \frac{\partial f\left(x^{*}\right)}{\partial x_{i}}\left(x_{i}-x_{i}^{*}\right) \geq 0, \quad \forall x_{i} \geq 0, i=1, \ldots, n
$$

Let us fix $i$. By letting $x_{j}=x_{j}^{*}$ for $j \neq i$ and $x_{i}=x_{i}^{*}+1$ in Eq. (1.2), we obtain

$$
\frac{\partial f\left(x^{*}\right)}{\partial x_{i}} \geq 0, \quad \forall i
$$

Furthermore, if we have $x_{i}^{*}>0$, by letting $x_{j}=x_{j}^{*}$ for $j \neq i$ and $x_{i}=\frac{1}{2} x_{i}^{*}$ in Eq. (1.2), we obtain $\partial f\left(x^{*}\right) / \partial x_{i} \leq 0$, which when combined with Eq. (1.3) yields

$$
\frac{\partial f\left(x^{*}\right)}{\partial x_{i}}=0, \quad \text { if } x_{i}^{*}>0
$$

---

![[ch_2_p6_img4.jpeg]]

Figure 2.1.4. Illustration of the optimality condition for an orthant constraint. At a minimum $x^{*}$, all the partial derivatives $\partial f\left(x^{*}\right) / \partial x_{i}$ are nonnegative, and they are zero for the inactive constraint indices, i.e. the indices $i$ with $x_{i}^{*}>0$. If all constraints are inactive, the optimality condition reduces to the unconstrained optimization condition $\nabla f\left(x^{*}\right)=0$.

The necessary conditions of Eqs. (1.3) and (1.4) are illustrated in Fig. 2.1.4.
Consider also the case where the constraints are upper and lower bounds on the variables, that is,

$$
X=\left\{x \mid \alpha_{i} \leq x_{i} \leq \beta_{i}, i=1, \ldots, n\right\}
$$

where $\alpha_{i}$ and $\beta_{i}$ are given scalars. We leave it as an exercise for the reader to verify that if $x^{*}$ is a local minimum, then

$$
\begin{aligned}
& \frac{\partial f\left(x^{*}\right)}{\partial x_{i}} \geq 0, \quad \text { if } x_{i}^{*}=\alpha_{i} \\
& \frac{\partial f\left(x^{*}\right)}{\partial x_{i}} \leq 0, \quad \text { if } x_{i}^{*}=\beta_{i} \\
& \frac{\partial f\left(x^{*}\right)}{\partial x_{i}}=0, \quad \text { if } \alpha_{i}<x_{i}^{*}<\beta_{i}
\end{aligned}
$$

Finally, suppose that $f$ is convex. Then, it can be seen that the conditions (1.3) and (1.4) imply the condition (1.2), which, by Prop. 2.1.2(b), is sufficient for optimality of $x^{*}$ in the case of nonnegativity constraints. Similarly, the conditions (1.6)-(1.8) are sufficient for optimality in the case of a convex $f$, and the upper and lower bound constraints of Eq. (1.5).

# Example 1.2 (Optimization Over a Simplex) 

Consider the case where the constraint set is a simplex

$$
X=\left\{x \mid x \geq 0, \sum_{i=1}^{n} x_{i}=r\right\}
$$

---

where $r>0$ is a given scalar. Then the necessary condition (1.1) for $x^{*}=$ $\left(x_{1}^{*}, \ldots, x_{n}^{*}\right)$ to be a local minimum is

$$
\sum_{i=1}^{n} \frac{\partial f\left(x^{*}\right)}{\partial x_{i}}\left(x_{i}-x_{i}^{*}\right) \geq 0, \quad \forall x_{i} \geq 0 \text { with } \sum_{i=1}^{n} x_{i}=r
$$

Fix an index $i$ for which $x_{i}^{*}>0$ and let $j$ be any other index. By using the feasible vector $x=\left(x_{1}, \ldots, x_{n}\right)$ with $x_{i}=0, x_{j}=x_{j}^{*}+x_{i}^{*}$, and $x_{m}=x_{m}^{*}$ for all $m \neq i, j$ in Eq. (1.10), we obtain

$$
\left(\frac{\partial f\left(x^{*}\right)}{\partial x_{j}}-\frac{\partial f\left(x^{*}\right)}{\partial x_{i}}\right) x_{i}^{*} \geq 0
$$

or equivalently

$$
x_{i}^{*}>0 \quad \Longrightarrow \quad \frac{\partial f\left(x^{*}\right)}{\partial x_{i}} \leq \frac{\partial f\left(x^{*}\right)}{\partial x_{j}}, \quad \forall j
$$

Thus, all coordinates which are positive at the optimum must have minimal (and equal) partial cost derivatives.

Assuming that $f$ is convex, we can show that Eq. (1.11) is also sufficient for global optimality of $x^{*}$. Indeed, suppose that $x^{*}$ is feasible and satisfies Eq. (1.11), and let

$$
\Delta=\min _{i=1, \ldots, n} \frac{\partial f\left(x^{*}\right)}{\partial x_{i}}
$$

For every $x \in X$, we have $\sum_{i=1}^{n}\left(x_{i}-x_{i}^{*}\right)=0$, so that

$$
0=\sum_{i=1}^{n} \Delta\left(x_{i}-x_{i}^{*}\right) \leq \sum_{\left\{i \mid x_{i}>x_{i}^{*}\right\}} \frac{\partial f\left(x^{*}\right)}{\partial x_{i}}\left(x_{i}-x_{i}^{*}\right)+\sum_{\left\{i \mid x_{i}<x_{i}^{*}\right\}} \Delta\left(x_{i}-x_{i}^{*}\right)
$$

If $i$ is such that $x_{i}<x_{i}^{*}$, we must have $x_{i}^{*}>0$ and, by condition (1.11), $\Delta=\partial f\left(x^{*}\right) / \partial x_{i}$. Thus $\Delta$ can be replaced by $\partial f\left(x^{*}\right) / \partial x_{i}$ in the right-hand side of the preceding inequality, thereby yielding Eq. (1.10), which by Prop. 2.1.2(b), implies that $x^{*}$ is optimal. Thus, when $f$ is convex, condition (1.11) is necessary and sufficient for optimality of $x^{*}$.

There is a straightforward generalization of this example to the case where $X$ is a Cartesian product of several simplices. Then, there is a separate condition of the form (1.11) for each simplex, that is, the condition

$$
\frac{\partial f\left(x^{*}\right)}{\partial x_{i}} \leq \frac{\partial f\left(x^{*}\right)}{\partial x_{j}}
$$

holds for all $i$ with $x_{i}^{*}>0$ and all $j$ for which $x_{j}$ is constrained by the same simplex as $x_{i}$.

---

# Example 1.3 (Optimal Routing in a Communication Network) 

We are given a directed graph, which is viewed as a model of a data communication network. We are also given a set $W$ of ordered node pairs $w=(i, j)$. The nodes $i$ and $j$ are referred to as the origin and the destination of $w$, respectively, and $w$ is referred to as an OD pair. For each $w$, we are given a scalar $r_{w}$ referred to as the input traffic of $w$. In the context of routing of data in a communication network, $r_{w}$ (measured in data units/second) is the arrival rate of traffic entering and exiting the network at the origin and the destination of $w$, respectively. The routing objective is to divide each $r_{w}$ among the many paths from origin to destination in a way that the resulting total arc flow pattern minimizes a suitable cost function. We denote:
$P_{w}$ : A given set of paths that start at the origin and end at the destination of $w$. All arcs on each of these paths are oriented in the direction from the origin to the destination.
$x_{p}$ : The portion of $r_{w}$ assigned to path $p$, also called the flow of path $p$.
The collection of all path flows $\left\{x_{p} \mid w \in W, p \in P_{w}\right\}$ must satisfy the constraints

$$
\begin{gathered}
\sum_{p \in P_{w}} x_{p}=r_{w}, \quad \forall w \in W \\
x_{p} \geq 0, \quad \forall p \in P_{w}, w \in W
\end{gathered}
$$

as shown in Fig. 2.1.5. The total flow $F_{i j}$ of arc $(i, j)$ is the sum of all path flows traversing the arc:

$$
F_{i j}=\sum_{\substack{\text { all paths } p \\
\text { containing }(i, j)}} x_{p}
$$

Consider a cost function of the form

$$
\sum_{(i, j)} D_{i j}\left(F_{i j}\right)
$$

The problem is to find a set of path flows $\left\{x_{p}\right\}$ that minimize this cost function subject to the constraints of Eqs. (1.12)-(1.14). We assume that $D_{i j}$ is a convex and continuously differentiable function of $F_{i j}$ with first derivative denoted by $D_{i j}^{\prime}$. In data routing applications, the form of $D_{i j}$ is often based on a queueing model of average delay (see [BeG92]).

The preceding problem is known as a multicommodity network flow problem. The terminology reflects the fact that the arc flows consist of several different commodities: in the present example, the different commodities are the data of the distinct OD pairs.

By expressing the total flows $F_{i j}$ in terms of the path flows in the cost function (1.15) [using Eq. (1.14)], the problem can be formulated in terms of the path flow variables $\left\{x_{p} \mid p \in P_{w}, w \in W\right\}$ as

$$
\begin{aligned}
& \text { minimize } D(x) \\
& \text { subject to } \sum_{p \in P_{w}} x_{p}=r_{w}, \quad \forall w \in W \\
& x_{p} \geq 0, \quad \forall p \in P_{w}, w \in W
\end{aligned}
$$

---

![[ch_2_p9_img5.jpeg]]

Figure 2.1.5. Constraints for the path flows of an OD pair $w$. The path flows $x_{p}$ of the paths $p \in P_{w}$ should be nonnegative and add up to the given traffic input $r_{w}$ of the OD pair.
where

$$
D(x)=\sum_{(i, j)} D_{i j}\left(\sum_{\substack{\text { all paths } p \\ \text { containing }(i, j)}} x_{p}\right)
$$

and $x$ is the vector of path flows $x_{p}$.
This problem, viewed as a problem in the variables $\left\{x_{p}\right\}$, has a convex and differentiable cost function and a constraint set that is a Cartesian product of simplices. Therefore, we can apply the condition [1.11, developed in the preceding example, separately, for the path flows of each OD pair. In particular, we will show that optimal routing directs traffic exclusively along paths that are shortest with respect to are lengths that depend on the flows carried by the arcs. To this end, we first derive the partial derivative of $D$ with respect to $x_{p}$. It is given by

$$
\frac{\partial D(x)}{\partial x_{p}}=\sum_{\substack{\text { all arcs }(i, j) \\ \text { on path } p}} D_{i j}^{\prime}\left(F_{i j}\right)
$$

where the first derivatives $D_{i j}^{\prime}$ are evaluated at the are flows $F_{i j}$ corresponding to $x$. From this equation, it is seen that $\partial D / \partial x_{p}$ is the length of path $p$ when the length of each arc (i, j) is taken to be the first derivative $D_{i j}^{\prime}$ evaluated at $F_{i j}$. Consequently, in what follows, $\partial D / \partial x_{p}$ is called the first derivative length of path $p$.

By using the necessary and sufficient optimality condition (1.11) of the preceding example, we obtain for all $w \in W$ and $p \in P_{w}$

$$
x_{p}^{*}>0 \quad \Longrightarrow \quad \frac{\partial D\left(x^{*}\right)}{\partial x_{p}} \leq \frac{\partial D\left(x^{*}\right)}{\partial x_{p^{\prime}}}, \quad \forall p^{\prime} \in P_{w}
$$

In words, the above condition says that a set of path flows is optimal if and only if path flow is positive only on paths with a minimum first derivative length. The condition also implies that at an optimum, the paths along which the input flow $r_{w}$ of OD pair $w$ is split must have equal length (and less than or equal length to that of all other paths of $w$ ).

---

# Example 1.4 (Traffic Assignment) 

We are given a directed graph, which is viewed as a model of a transportation network. The arcs of the graph represent transportation links such as highways, rail lines, etc. The nodes of the graph represent junction points where traffic can exit from one transportation link and enter another. Similar to the preceding example, we are given a set $W$ of OD pairs. For OD pair $w=(i, j)$, there is a known input $r_{w}>0$ representing traffic entering the network at the origin node $i$ of $w$ and exiting the network at the destination node $j$ of $w$. The input $r_{w}$ is to be divided among a set $P_{w}$ of paths starting at the origin node of $w$ and ending at the destination node of $w$. Let $x_{p}$ denote the portion of $r_{w}$ carried by path $p$ and let $x$ be the vector having as coordinates all the path flows $x_{p}, p \in P_{w}, w \in W$. Thus, the constraint set is the Cartesian product of simplices defined by Eqs. (1.12) and (1.13) of the preceding example.

For each arc $(i, j)$, we are given a function $T_{i j}\left(F_{i j}\right)$ of the total traffic $F_{i j}$ carried by arc $(i, j)$ [cf. Eq. (1.14)]. This function models the time required for traffic to travel from the start node to the end node of the arc $(i, j)$. An important problem is to find a path flow vector $x^{*} \in X$ that consists of path flows that are positive only on paths of minimum travel time. That is, for all $w \in W$ and paths $p \in P_{w}$, we require that

$$
x_{p}^{*}>0 \quad \Longrightarrow \quad t_{p}\left(x^{*}\right) \leq t_{p^{\prime}}\left(x^{*}\right), \quad \forall p^{\prime} \in P_{w}, \forall w \in W
$$

where $t_{p}(x)$, the travel time of path $p$, is defined as the sum of the travel times of the arcs of the path,

$$
t_{p}(x)=\sum_{\substack{\text { all arcs }(i, j) \\ \text { on path } p}} T_{i j}\left(F_{i j}\right), \quad \forall p \in P_{w}, \forall w \in W
$$

The preceding problem draws its validity from a hypothesis, called the user-optimization principle, which asserts that traffic equilibrium is established when each user of the network chooses, among all available paths, a path requiring minimum travel time. Thus, assuming that the user-optimization principle holds, a path flow vector $r^{*}$ that solves the problem also models accurately the distribution of traffic through the network, and can be used for traffic projections when planning modifications to the network.

We now observe that the minimum travel time hypothesis (1.19) is identical with the optimality condition (1.18) if we identify the arc travel time $T_{i j}\left(F_{i j}\right)$ with the cost derivative $D_{i j}^{\prime}\left(F_{i j}\right)$ [cf. Eq. (1.17)]. It follows that we can solve the transportation problem by converting it to the optimal routing problem of the preceding example with the identification

$$
D_{i j}\left(F_{i j}\right)=\int_{0}^{F_{i j}} T_{i j}(\xi) d \xi
$$

If we assume that $T_{i j}$ is continuous and monotonically nondecreasing, as is natural in a transportation context, it is straightforward to show that the function $D_{i j}$ as defined by Eq. (1.20) is convex with derivative $D_{i j}^{\prime}$ equal to $T_{i j}$, so a minimum first derivative length path is a path of minimum travel time.

---

# Projection on a Convex Set 

Let $z$ be a fixed vector in $\Re^{n}$ and consider the problem of finding a vector $x^{*}$ of the closed convex set $X$, which is at a minimum distance from $z$; that is,

$$
\begin{array}{ll}
\operatorname{minimize} & f(x)=\|z-x\|^{2} \\
\text { subject to } & x \in X
\end{array}
$$

We call this the problem of projection of $z$ on $X$. This problem has important applications in mathematical analysis and plays also an important role in nonlinear programming algorithms. The main facts are summarized in the following proposition.

## Proposition 2.1.3: (Projection Theorem)

(a) For every $z \in \Re^{n}$, there exists a unique $x^{*} \in X$ that minimizes $\|z-x\|$ over all $x \in X$. This vector is called the projection of $z$ on $X$ and is denoted by $[z]^{+}$.
(b) Given some $z \in \Re^{n}$, a vector $x^{*} \in X$ is equal to the projection $[z]^{+}$if and only if

$$
\left(z-x^{*}\right)^{\prime}\left(x-x^{*}\right) \leq 0, \quad \forall x \in X
$$

(See Fig. 2.1.6.)
(c) The mapping $f: \Re^{n} \mapsto X$ defined by $f(x)=[x]^{+}$is continuous and nonexpansive, that is,

$$
\left\|[x]^+-[y]^{+}\right\| \leq\|x-y\|, \quad \forall x, y \in \Re^{n}
$$

(d) In the case where $X$ is a [[Subspaces|subspace]], a vector $x^{*} \in X$ is equal to the projection $[z]^{+}$if and only if $z-x^{*}$ is orthogonal to $X$, that is,

$$
\left(z-x^{*}\right)^{\prime} x=0, \quad \forall x \in X
$$

Proof: Parts (a) and (c) are proved in [[Mathematics/Optimization/Nonlinear Optimization/Textbook/Appendix#The Projection Theorem|Prop. B. 11]] of Appendix B. Part (b) is also proved in Prop. B. 11 and is equivalent to the necessary and sufficient condition of [[Chapter 2#Proposition 2.1.2 (Optimality Condition)|Prop. 2.1.2]], specialized to the projection problem. To prove part (d), note that since  is a subspace, the vectors $x^{*}+x$ and $x^{*}-x$ belong to  for all $x \in X$. Applying the condition of part (b) with these vectors replacing $x$, we obtain $\left(z-x^{*}\right)^{\prime} x \leq 0$ and $\left(z-x^{*}\right)^{\prime}(-x) \leq 0$, which imply the desired result. Q.E.D.

---

![[Mathematics/Optimization/Nonlinear Optimization/Textbook/Images/ch_2_p12_img6.jpeg]]

Figure 2.1.6. Illustration of the necessary and sufficient condition for $x^{*}$ to be the unique projection of $z$ on the closed convex set $X$. The angle between $z-x^{*}$ and $x-x^{*}$ should be greater or equal to 90 degrees for all $x \in X$.

# Example 1.5 (Quadratic Programming with Equality Constraints) 

Consider the quadratic programming problem

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2}\|x\|^{2}+c^{\prime} x \\
& \text { subject to } A x=0
\end{aligned}
$$

where $c$ is a given vector in $\Re^{n}$ and $A$ is an $m \times n$ matrix of rank $m$.
By adding the constant term $\frac{1}{2}\|c\|^{2}$ to the cost function, we can equivalently write this problem as

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2}\|c+x\|^{2} \\
& \text { subject to } A x=0
\end{aligned}
$$

which is the problem of projecting the vector $-c$ on the subspace $X=\{x \mid$ $A x=0\}$. By Prop. 2.1.3(d), a vector $x^{*}$ such that $A x^{*}=0$ is the unique projection if and only if

$$
\left(c+x^{*}\right)^{\prime} x=0, \quad \forall x \text { with } A x=0
$$

It can be seen that the vector

$$
x^{*}=-\left(I-A^{\prime}\left(A A^{\prime}\right)^{-1} A\right) c
$$

satisfies this condition and is thus the unique solution of the quadratic programming problem (1.23). [Note here that the matrix $A A^{\prime}$ is invertible because $A$ has rank $m$ (Prop. A. 20 in Appendix A).]

Note that $x^{*}$ is zero if and only if $c$ is orthogonal to the subspace $X=\{x \mid A x=0\}$, or equivalently from Eq. (1.24),

$$
c=A^{\prime}\left(A A^{\prime}\right)^{-1} A c
$$

---

Thus, in this case, $c$ can be expressed in the form $A^{\prime} \mu$, where $\mu \in \Re^{m}$, that is, as a linear combination of the rows of $A$. Conversely, it is seen that every vector of the form $A^{\prime} \mu$, with $\mu \in \Re^{m}$, is orthogonal to the subspace $X$.

The conclusion is that the orthogonal complement $X^{\perp}$ of $X$ (the set of all vectors orthogonal to all $x$ with $A x=0$ ), is the subspace $\{x \mid x=$ $\left.A^{\prime} \mu, \mu \in \Re^{m}\right\}$. Our proof of this fact assumed that $A$ has rank $m$ but can be easily modified for the case where the rows of $A$ are linearly dependent by eliminating a sufficient number of dependent rows of $A$. Also a more general version of the preceding analysis, involving cones rather than subspaces, is given in Prop. B. 16 of Appendix B.

Consider now the more general quadratic program

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2}(x-\bar{x})^{\prime} H(x-\bar{x})+c^{\prime}(x-\bar{x}) \\
& \text { subject to } A x=b
\end{aligned}
$$

where $c$ and $A$ are as before, $H$ is a positive definite symmetric matrix, $b$ is a given vector in $\Re^{m}$, and $\bar{x}$ is a given vector in $\Re^{n}$, which is feasible, that is, satisfies $A \bar{x}=b$. By introducing the vector $y=H^{1 / 2}(x-\bar{x})$, we can write this problem as

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2}\|y\|^{2}+\left(H^{-1 / 2} c\right)^{\prime} y \\
& \text { subject to } A H^{-1 / 2} y=0
\end{aligned}
$$

Using Eq. (1.24) we see that the solution of this problem is

$$
y^{*}=-\left(I-H^{-1 / 2} A^{\prime}\left(A H^{-1} A^{\prime}\right)^{-1} A H^{-1 / 2}\right) H^{-1 / 2} c
$$

and by passing to the $x$-coordinate system through the transformation $x^{*}-$ $\bar{x}=H^{-1 / 2} y^{*}$, we obtain the optimal solution

$$
x^{*}=\bar{x}-H^{-1}\left(c-A^{\prime} \lambda\right)
$$

where the vector $\lambda$ is given by

$$
\lambda=\left(A H^{-1} A^{\prime}\right)^{-1} A H^{-1} c
$$

The preceding program (1.25) contains as a special case the positive definite quadratic program

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2} x^{\prime} H x+c^{\prime} x \\
& \text { subject to } A x=b
\end{aligned}
$$

This special case is obtained when $\bar{x}$ is given by

$$
\bar{x}=H^{-1} A^{\prime}\left(A H^{-1} A^{\prime}\right)^{-1} b
$$

Indeed $\bar{x}$ as given above satisfies $A \bar{x}=b$ as required, and for all $x$ with $A x=b$, we have

$$
x^{\prime} H \bar{x}=x^{\prime} A^{\prime}\left(A H^{-1} A^{\prime}\right)^{-1} b=b^{\prime}\left(A H^{-1} A^{\prime}\right)^{-1} b
$$

---

which implies that for all $x$ with $A x=b$
$\frac{1}{2}(x-\bar{x})^{\prime} H(x-\bar{x})+c^{\prime}(x-\bar{x})=\frac{1}{2} x^{\prime} H x+c^{\prime} x+\left(\frac{1}{2} \bar{x}^{\prime} H \bar{x}-c^{\prime} \bar{x}-b^{\prime}\left(A H^{-1} A^{\prime}\right)^{-1} b\right)$.
The last term in parentheses on the right-hand side above is constant, thus establishing that the programs (1.25) and (1.28) have the same optimal solution when $\bar{x}$ is given by Eq. (1.29). By combining Eqs. (1.26) and 1.29), we obtain the optimal solution of program (1.28):

$$
x^{*}=-H^{-1}\left(c-A^{\prime} \lambda-A^{\prime}\left(A H^{-1} A^{\prime}\right)^{-1} b\right)
$$

where $\lambda$ is given by Eq. (1.27).

# E X E R C ISES 

### 1.1 (Heron's Problem)

In three-dimensional space, consider a two-dimensional plane and two points $z_{1}$ and $z_{2}$ lying outside the plane. Use the optimality condition of Prop. 2.1.2 to characterize the vector $x^{*}$, which minimizes $\left\|z_{1}-x\right\|+\left\|z_{2}-x\right\|$ over all $x$ in the plane.

### 1.2 (Euclid's Problem)

Among all parallelograms that are contained in a given triangle, find the ones that have maximal area. Hint: First argue that an optimal solution can be found among the parallelograms that share one angle with the triangle.

### 1.3 (Archimedes' Problem)

A spherical segment is the intersection of a sphere with one of the halfspaces corresponding to a 2-dimensional plane that meets the sphere. The area of a spherical segment is the area of the portion of its surface that is also part of the surface of the sphere. Note that if $r$ is the radius of the sphere and $h$ is the height of a spherical segment, the volume enclosed by the segment is $\pi h^{2}(r-$ $h / 3)$ and its spherical area is $2 \pi r h$. Show that among all spherical segments with the same spherical area, the one that encloses the largest volume is a hemisphere.

---

# 1.4 (Kepler's Planimetric Problem) 

Among all rectangles contained in a given circle show that the one that has maximal area is a square.

### 1.5 (Tartaglia's Problem)

Divide the number 8 into two nonnegative parts $x$ and $y$ so as to maximize $x y(x-y)$. Answer: $x=4+4 / \sqrt{3}, y=4-4 / \sqrt{3}$.

## 1.6

Consider the problem

$$
\begin{aligned}
& \text { maximize } x_{1}^{a_{1}} x_{2}^{a_{2}} \cdots x_{n}^{a_{n}} \\
& \text { subject to } \sum_{i=1}^{n} x_{i}=1, \quad x_{i} \geq 0, \quad i=1, \ldots, n
\end{aligned}
$$

where $a_{i}$ are given positive scalars. Find a global maximum and show that it is unique.

## 1.7

Let $Q$ be a positive definite symmetric matrix. State and prove a generalization of the projection theorem that involves the cost function $(z-x)^{\prime} Q(z-x)$ in place of $\|z-x\|^{2}$. Hint: Use a transformation of variables.

## 1.8

Consider a triangle on a plane and a point $z$ outside the plane whose projection on the plane is $\hat{z}$. Characterize the projection of $z$ on the triangle in terms of the position of $\hat{z}$ relative to the triangle.

## 1.9

Verify the necessary optimality conditions of Eqs. (1.6)-(1.8). Show that, assuming $f$ is convex, they are also sufficient for optimality of $x^{*}$.

### 1.10 (Second Order Necessary Optimality Condition)

Show that if $x^{*}$ is a local minimum of the twice continuously differentiable function $f: \Re^{n} \mapsto \Re$ over the convex set $X$, then

$$
\left(x-x^{*}\right)^{\prime} \nabla^{2} f\left(x^{*}\right)\left(x-x^{*}\right) \geq 0
$$

for all $x \in X$ such that $\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right)=0$.

---

# 1.11 (Second Order Sufficiency Conditions) 

Show that a vector $x^{*} \in X$ is a local minimum of a twice continuously differentiable function $f: \Re^{n} \mapsto \Re$ over the convex set $X$ if

$$
\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right) \geq 0, \quad \forall x \in X
$$

and one of the following two conditions holds:
(1) $X$ is polyhedral and we have

$$
\left(x-x^{*}\right)^{\prime} \nabla^{2} f\left(x^{*}\right)\left(x-x^{*}\right)>0
$$

for all $x \in X$ satisfying $x \neq x^{*}$ and $\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right)=0$.
(2) For some $\gamma>0$, we have

$$
\left(x-x^{*}\right)^{\prime} \nabla^{2} f\left(x^{*}\right)\left(x-x^{*}\right) \geq \gamma\left\|x-x^{*}\right\|^{2} . \quad \forall x \in X
$$

### 1.12 (Projection on a Simplex)

(a) Develop an algorithm to find the projection of a vector $z$ on a simplex. This algorithm should be almost as simple as a closed form solution.
(b) Modify the algorithm of part (a) so that it finds the minimum of the cost function

$$
f(x)=\sum_{i=1}^{n}\left(\alpha_{i} x_{i}+\frac{1}{2} \beta_{i} x_{i}^{2}\right)
$$

over a simplex, where $\alpha_{i}$ and $\beta_{i}$ are given scalars with $\beta_{i}>0$. Consider also the case where $\beta_{i}$ may be zero for some indices $i$.

### 1.13 (Optimal Control - Minimum Principle)

Consider the problem of finding sequences $u=\left(u_{0}, u_{1}, \ldots, u_{N-1}\right)$ and $x=$ $\left(x_{1}, x_{2}, \ldots, x_{N}\right)$ that minimize

$$
g_{N}\left(x_{N}\right)+\sum_{i=0}^{N-1} g_{i}\left(x_{i}, u_{i}\right)
$$

subject to the constraints

$$
\begin{gathered}
x_{i+1}=f_{i}\left(x_{i}, u_{i}\right), \quad i=0, \ldots, N-1, \quad x_{0}: \text { given, } \\
u_{i} \in U_{i} \subset \Re^{m}, \quad i=0, \ldots, N-1
\end{gathered}
$$

---

All functions $g_{i}$ and $f_{i}$ are assumed continuously differentiable, and the sets $U_{i}$ are assumed convex. Show that if $u^{*}$ and $x^{*}$ are optimal, then

$$
\nabla_{u_{i}} H_{i}\left(x_{i}^{*}, u_{i}^{*}, p_{i+1}^{*}\right)^{\prime}\left(u_{i}-u_{i}^{*}\right) \geq 0, \quad \forall u_{i} \in U_{i}, i=0, \ldots, N-1
$$

where

$$
H_{i}\left(x_{i}, u_{i}, p_{i+1}\right)=g_{i}\left(x_{i}, u_{i}\right)+p_{i+1}^{\prime} f_{i}\left(x_{i}, u_{i}\right)
$$

is the Hamiltonian function, and the vectors $p_{1}^{*}, \ldots, p_{N}^{*}$ are obtained from the adjoint equation

$$
p_{i}^{*}=\nabla_{x_{i}} H_{i}\left(x_{i}^{*}, u_{i}^{*}, p_{i+1}^{*}\right), \quad i=1, \ldots, N-1
$$

with the terminal condition

$$
p_{N}^{*}=\nabla g_{N}\left(x_{N}^{*}\right)
$$

If, in addition, the Hamiltonian $H_{i}$ is a convex function of $u_{i}$ for any fixed $x_{i}$ and $p_{i+1}$, we have

$$
u_{i}^{*}=\arg \min _{u_{i} \in U_{i}} H_{i}\left(x_{i}^{*}, u_{i}, p_{i+1}^{*}\right), \quad \forall i=0, \ldots, N-1
$$

(The last condition is known as the minumum principle.) Hint: Apply Prop. 2.1.2 and use the results of Section 1.9.

# 1.14 

A farmer annually producing $x_{i}$ units of a certain crop stores $\left(1-u_{i}\right) x_{i}$ units of his production, where $0 \leq u_{i} \leq 1$, and invests the remaining $u_{i} x_{i}$ units, thus increasing the next year's production to a level $x_{i+1}$ given by

$$
x_{i+1}=x_{i}+u u_{i} x_{i}, \quad i=0,1, \ldots, N-1
$$

where $u$ is a given positive scalar. The problem is to find the optimal investment sequence $u_{0}, \ldots, u_{N-1}$ that maximizes the total product stored over $N$ years

$$
x_{N}+\sum_{i=0}^{N-1}\left(1-u_{i}\right) x_{i}
$$

Show that one optimal sequence is given by:
(1) If $w>1, u_{0}^{*}=\cdots=u_{N-1}^{*}=1$.
(2) If $0<w<1 / N, u_{0}^{*}=\cdots=u_{N-1}^{*}=0$.
(3) If $1 / N \leq w \leq 1$,

$$
\begin{aligned}
& u_{0}^{*}=\cdots=u_{N-\bar{i}-1}=1 \\
& u_{N-\bar{i}}^{*}=\cdots=u_{N-1}^{*}=0
\end{aligned}
$$

where $\bar{i}$ is such that $1 /(\bar{i}+1)<w \leq 1 / \bar{i}$.

---

# 1.15 

Let $x_{i}$ denote the number of educators in a certain country at time $i$ and let $y_{i}$ denote the number of research scientists at time $i$. New scientists (potential educators or research scientists) are produced during the $i$ th period by educators at a rate $\gamma$ per educator, while educators and research scientists leave the field due to death, retirement, and transfer at a rate $\delta$, where $\gamma$ and $\delta$ are given scalars with $0<\gamma$ and $0<\delta<1$. By means of incentives a science policy maker can determine the proportion $u_{i}$ of new scientists produced at time $i$ who become educators. Thus the number of research scientists and educators evolves according to the equations

$$
\begin{gathered}
x_{i+1}=(1-\delta) x_{i}+u_{i} \gamma x_{i} \\
y_{i+1}=(1-\delta) y_{i}+\left(1-u_{i}\right) \gamma x_{i}
\end{gathered}
$$

The initial numbers $x_{0}, y_{0}$ are known, and the $u_{i}$ are constrained by

$$
0<\alpha \leq u_{i} \leq \beta<1, \quad i=0,1, \ldots, N-1
$$

where the scalars $\alpha$ and $\beta$ are given. Find a sequence $\left\{u_{0}^{*}, \ldots, u_{N-1}^{*}\right\}$ that maximizes the final number $y_{N}$ of research scientists.

### 1.16 (Fractional Programming)

Consider the problem

$$
\begin{aligned}
& \operatorname{minimize} \frac{f(x)}{g(x)} \\
& \text { subject to } x \in X
\end{aligned}
$$

where $f: \Re^{n} \mapsto \Re$ and $g: \Re^{n} \mapsto \Re$ are given functions and $X$ is a given subset such that $g(x)>0$ for all $x \in X$. For $\lambda \in \Re$, define

$$
Q(\lambda)=\min _{x \in X}\{f(x)-\lambda g(x)\}
$$

and suppose that a scalar $\lambda^{*}$ and a vector $x^{*} \in X$ satisfy $Q\left(\lambda^{*}\right)=0$ and

$$
x^{*}=\arg \min _{x \in X}\left\{f(x)-\lambda^{*} g(x)\right\}
$$

Show that $x^{*}$ is an optimal solution of the original problem.

---

# 2.2 FEASIBLE DIRECTIONS AND THE CONDITIONAL GRADIENT METHOD 

We now turn to computational methods for solving the constrained problem of the preceding section:

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } x \in X
\end{aligned}
$$

We will see that there is a great variety of algorithms for this problem. However, in this chapter we will restrict ourselves to a limited class of methods that have the following characteristics:
(a) They do not rely on any structure of the constraint set other than its convexity.
(b) They generate sequences of feasible points $\left\{x^{k}\right\}$ by searching along descent directions. In this sense, they can be viewed as constrained versions of the unconstrained descent algorithms of the previous chapter.

Most of the algorithms in this chapter belong to the class of the, so called, feasible direction methods, which we proceed to discuss.

### 2.2.1 Descent Directions and Stepsize Rules

Given a feasible vector $x^{k}$, a feasible direction at $x^{k}$ is a vector $d^{k} \neq 0$ such that $x^{k}+\alpha d^{k}$ is feasible for all $\alpha>0$ that are sufficiently small. Figure 2.2.1 illustrates the set of feasible directions at a point.
![[ch_2_p19_img7.jpeg]]

Figure 2.2.1. Feasible directions $d$ at a feasible point $x$. By definition, $d$ is a feasible direction if changing $x$ by a small amount in the direction $d$ maintains feasibility.

---

A feasible direction method starts with a feasible vector $x^{0}$ and generates a sequence of feasible vectors $\left\{x^{k}\right\}$ according to

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}
$$

where, if $x^{k}$ is not stationary, $d^{k}$ is a feasible direction such that

$$
\nabla f\left(x^{k}\right)^{\prime} d^{k}<0
$$

and the stepsize $\alpha^{k}$ is chosen to be positive and such that

$$
x^{k}+\alpha^{k} d^{k} \in X
$$

If $x^{k}$ is stationary, the method stops, that is, $x^{k+1}=x^{k}$ (equivalently, we choose $d^{k}=0$ ). We will primarily concentrate on feasible direction methods that are also descent algorithms; that is, the stepsize $\alpha^{k}$ is selected so that

$$
f\left(x^{k}+\alpha^{k} d^{k}\right)<f\left(x^{k}\right), \quad \forall k
$$

Figure 2.2.2 illustrates a feasible direction method.
![[ch_2_p20_img8.jpeg]]

Figure 2.2.2. Sample path of a feasible direction method. At each iteration, a feasible direction of descent $d^{k}$ at the current point $x^{k}$ is chosen and a new feasible point $x^{k+1}=x^{k}+\alpha^{k} d^{k}$ is obtained along $d^{k}$.

In our case where $X$ is convex, there is an alternative and equivalent characterization of feasible direction methods. In particular, it can be seen that the feasible directions at $x^{k}$ are the vectors of the form

$$
d^{k}=\gamma\left(\bar{x}^{k}-x^{k}\right), \quad \gamma>0
$$

where $\bar{x}^{k}$ is some feasible vector different from $x^{k}$. Thus, a feasible direction method can be written in the form

$$
x^{k+1}=x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right)
$$

---

where

$$
\alpha^{k} \in(0,1]
$$

and if $x^{k}$ is nonstationary,

$$
\bar{x}^{k} \in X, \quad \nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right)<0
$$

Note that because the constraint set $X$ is convex, we have $x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right) \in$ $X$ for all $\alpha^{k} \in[0,1]$ when $x^{k} \in X$, so that the generated sequence of iterates $\left\{x^{k}\right\}$ is feasible. Furthermore, if $x^{k}$ is nonstationary, there always exists a feasible direction $\bar{x}^{k}-x^{k}$ with the descent property (2.7), since otherwise we would have $\nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right) \geq 0$ for all $x \in X$, implying that $x^{k}$ is stationary.

# Stepsize Selection in Feasible Direction Methods 

Most of the rules for choosing the stepsize $\alpha^{k}$ in gradient methods apply also to feasible direction methods. Some of the most popular ones are:

## Limited Minimization Rule

Here, $\alpha^{k}$ is chosen so that

$$
f\left(x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right)\right)=\min _{\alpha \in[0,1]} f\left(x^{k}+\alpha\left(\bar{x}^{k}-x^{k}\right)\right)
$$

Note that there is no loss of generality in limiting the stepsize to the interval $[0,1]$, since different stepsize ranges can in effect be used by redefining the feasible vector $\bar{x}^{k}$.

## Armijo Rule

Here, fixed scalars $\beta$, and $\sigma>0$, with $\beta \in(0,1)$, and $\sigma \in(0,1)$ are chosen, and we set $\alpha^{k}=\beta^{m_{k}} s$, where $m_{k}$ is the first nonnegative integer $m$ for which

$$
f\left(x^{k}\right)-f\left(x^{k}+\beta^{m}\left(\bar{x}^{k}-x^{k}\right)\right) \geq-\sigma \beta^{m} \nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right)
$$

In other words, the stepsizes $1, \beta, \beta^{2}, \ldots$, are tried successively until the above inequality is satisfied for $m=m_{k}$.

## Constant Stepsize

Here, a fixed stepsize

$$
\alpha^{k}=1, \quad k=0,1, \ldots
$$

is used. This choice is not as simple or as restrictive as it may seem. In fact, in any feasible direction method, we can use a constant unity stepsize if we rescale or redefine appropriately the direction $\bar{x}^{k}-x^{k}$. In this case, the burden is placed in effect on the direction finding procedure to yield a direction that guarantees descent and convergence.

---

# Convergence Analysis of Feasible Direction Methods 

The convergence analysis of feasible direction methods for the case of the minimization rule and the Armijo rule is very similar to the one for gradient methods. As in Section 1.2, we say that the direction sequence $\left\{d^{k}\right\}$ is gradient related to $\left\{x^{k}\right\}$ if the following property can be shown:

For any subsequence $\left\{x^{k}\right\}_{k \in \mathcal{K}}$ that converges to a nonstationary point, the corresponding subsequence $\left\{d^{k}\right\}_{k \in \mathcal{K}}$ is bounded and satisfies

$$
\limsup _{k \rightarrow \infty, k \in \mathcal{K}} \nabla f\left(x^{k}\right)^{\prime} d^{k}<0
$$

A verbatim repetition of the proof of Prop. 1.2.1 of Section 1.2 shows the following result:

**Proposition 2.2.1**: (Stationarity of Limit Points for Feasible Direction Methods) Let $\left\{x^{k}\right\}$ be a sequence generated by the feasible direction method $x^{k+1}=x^{k}+\alpha^{k} d^{k}$. Assume that $\left\{d^{k}\right\}$ is gradient related and that $\alpha^{k}$ is chosen by the limited minimization rule or the Armijo rule. Then every limit point of $\left\{x^{k}\right\}$ is a stationary point. ^bb845e

There are also convergence results for feasible direction methods in the case of a constant stepsize, but they are best considered in the context of individual methods.

## Finding an Initial Feasible Point

To apply a feasible direction method, it is necessary to have an initial feasible point. Finding such a point may present difficulties if the constraint set is specified by nonlinear inequality constraints. If, however, the constraint set is a polyhedron (is specified by linear equality and inequality constraints), one can find an initial feasible point by solving the linear programming problem involving the same constraints and an arbitrary linear cost function.

An alternative method is based on introducing an additional artificial variable $y$, which is constrained to be nonnegative but is heavily penalized in the cost so that it is zero at the optimum. In particular, the problem

$$
\begin{aligned}
& \text { minimize } f(x) \\
& \text { subject to } a_{i}^{\prime} x=b_{i}, \quad x_{i} \geq 0, \quad i=1, \ldots, n
\end{aligned}
$$

---

is replaced by the problem
minimize $f(x)+c y$
subject to $a_{i}^{\prime} x+\left(b_{i}-\sum_{j=1}^{n} a_{i j}\right) y=b_{i}, \quad x_{i} \geq 0, \quad i=1, \ldots, n, \quad y \geq 0$,
where $c$ is a very large cost coefficient, $a_{i j}$ is the $j$ th coordinate of the vector $a_{i}$. The vector $(\bar{x}, \bar{y})$ with $\bar{x}_{i}=1, i=1, \ldots, n, \bar{y}=1$ is feasible for the modified problem. Note that this vector is interior to the bound constraints $x_{i} \geq 0, y \geq 0$; this is significant for interior point methods, which will be discussed in Section 2.6.

Similarly, the problem

$$
\begin{aligned}
& \text { minimize } f(x) \\
& \text { subject to } a_{j}^{\prime} x \leq b_{j}, \quad j=1, \ldots, r
\end{aligned}
$$

can be converted to the problem

$$
\begin{aligned}
& \text { minimize } f(x)+c y \\
& \text { subject to } a_{j}^{\prime} x-y \leq b_{j}, \quad j=1, \ldots, r, \quad y \geq 0
\end{aligned}
$$

where $c$ is a very large cost coefficient. For any vector $\bar{x}$ that is infeasible for the original problem, the vector $(\bar{x}, \bar{y})$ where

$$
\bar{y}=\max _{j=1, \ldots, r}\left\{a_{j}^{\prime} \bar{x}-b_{j}\right\}
$$

is feasible for the modified problem, while the vector $(\bar{x}, \bar{y}+1)$ is feasible, as well as interior to the inequality constraints of the modified problem.

# 2.2.2 The Conditional Gradient Method 

The most straightforward way to generate a feasible direction $\bar{x}^{k}-$ $x^{k}$ satisfying the descent condition $\nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right)<0$ is to solve the optimization problem

$$
\begin{aligned}
& \operatorname{minimize} \nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right) \\
& \text { subject to } x \in X
\end{aligned}
$$

and obtain $\bar{x}^{k}$ as the solution, that is,

$$
\bar{x}^{k}=\arg \min _{x \in X} \nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)
$$

Here we assume that $X$ is compact, so that the direction generating subproblem (2.11) has a solution. The corresponding feasible direction method

---

![[ch_2_p24_img9.jpeg]]

Figure 2.2.3. Finding the feasible direction $\bar{x}-x$ at a point $x$ in the conditional gradient method. $\bar{x}$ is a point of $X$ that lies furthest along the negative gradient direction $-\nabla f(x)$.

is the, so called, **conditional gradient method**, also known as the **FrankWolfe method**. The process to obtain $\bar{x}^{k}$ is illustrated in Fig. 2.2.3. In particular, $\bar{x}^{k}$ is the remotest point of $X$ along the negative gradient direction.

Naturally, in order for the method to make practical sense, the subproblem (2.11) must be much simpler than the original. This is typically the case when $f$ is nonlinear, and $X$ is specified by linear equality and inequality constraints. Then, the subproblem (2.11) is a linear program, which in many applications is very easy to solve. A principal example is when $X$ is a simplex or a Cartesian product of simplices, as in the communication and transportation examples of the previous section.

# Example 2.1 (Simplex Constraint) 

Let

$$
X=\left\{x \mid x \geq 0, \sum_{i=1}^{n} x_{i}=r\right\}
$$

where $r>0$ is a given scalar (cf. Example 1.2 of the previous section). Then, the subproblem (2.11) takes the form

$$
\begin{aligned}
& \operatorname{minimize} \sum_{i=1}^{n} \frac{\partial f\left(x^{k}\right)}{\partial x_{i}}\left(x_{i}-x_{i}^{k}\right) \\
& \text { subject to } \sum_{i=1}^{n} x_{i}=r, \quad x_{i} \geq 0
\end{aligned}
$$

and a solution $\bar{x}^{k}$ has all coordinates equal to zero except for a single coordinate, say $j$, which is equal to $r$; the $j$ th coordinate corresponds to a minimal partial derivative, that is,

$$
j=\arg \min _{i=1, \ldots, n} \frac{\partial f\left(x^{k}\right)}{\partial x_{i}}
$$

---

This choice, together with the path followed by the conditional gradient method is illustrated in Fig. 2.2.4.
![[ch_2_p25_img10.jpeg]]

Figure 2.2.4. Successive iterates of the conditional gradient method with the limited minimization stepsize rule, for the case of a simplex constraint.

# Convergence of the Conditional Gradient Method 

We will show that the direction sequence of the conditional gradient method is gradient related, so [[Chapter 2#^bb845e|Prop. 2.2.1]] applies. Indeed, suppose that $\left\{x^{k}\right\}_{k \in K}$ converges to a nonstationary point $\tilde{x}$. We must prove that

$$
\begin{gathered}
\limsup _{k \rightarrow \infty, k \in K}\left\|\bar{x}^{k}-x^{k}\right\|<\infty \\
\limsup _{k \rightarrow \infty, k \in K} \nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right)<0
\end{gathered}
$$

Equation (2.14) holds because $\bar{x}^{k} \in X, x^{k} \in X$, and the set $X$ is assumed compact. To show Eq. (2.15), we note that, by the definition of $\bar{x}^{k}$.

$$
\nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right) \leq \nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right), \quad \forall x \in X
$$

and by taking limit we obtain

$$
\limsup _{k \rightarrow \infty, k \in K} \nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right) \leq \nabla f(\tilde{x})^{\prime}(x-\tilde{x}), \quad \forall x \in X
$$

By taking the minimum over $x \in X$ and by using the nonstationarity of $\tilde{x}$, we have

$$
\limsup _{k \rightarrow \infty, k \in K} \nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right) \leq \min _{x \in X} \nabla f(\tilde{x})^{\prime}(x-\tilde{x})<0
$$

thereby proving Eq. (2.15).
We have thus shown that every limit point of the conditional gradient method with the limited minimization rule or the Armijo rule is stationary. An alternative line of convergence analysis is given in Exercise 2.4.

---

# Rate of Convergence of the Conditional Gradient Method 

Unfortunately, the asymptotic rate of convergence of the conditional gradient method is not very fast when the set $X$ is a polyhedron, i.e., it is specified by linear equality and inequality constraints. A partial explanation is that the vectors $\bar{x}^{k}$ used in the algorithm are typically extreme points (vertices) of $X$. Thus, the feasible direction used may tend to be orthogonal to the direction leading to the minimum (see Fig. 2.2.4). In fact one may construct simple examples where the error sequences $\left\{f\left(x^{k}\right)-f\left(x^{*}\right)\right\}$ and $\left\{\left\|x^{k}-x^{*}\right\|\right\}$ do not converge linearly (see Exercise 2.3). Thus, when $X$ is a polyhedron, the method is recommended only for problems where solution accuracy is not very important.

Somewhat peculiarly, the practical performance of the conditional gradient method tends to improve as the number of constraints becomes very large, even if these constraints are linear. An explanation for this is given in [Dun79], [DuS83], where it is shown that the convergence rate of the method is linear when the constraint set is not polyhedral but rather has a "positive curvature" property (for example it is a sphere). When there are many linear constraints, the constraint set tends to have very many closely spaced extreme points, and has this "positive curvature" property in an approximate sense.

## E XERCISES

## 2.1 (Computational Problem)

Consider the three-dimensional problem

$$
\begin{aligned}
& \operatorname{minimize} f(x)=\frac{1}{2}\left(x_{1}^{2}+x_{2}^{2}+0.1 x_{3}^{2}\right)+0.55 x_{3} \\
& \text { subject to } x_{1}+x_{2}+x_{3}=1, \quad 0 \leq x_{1}, 0 \leq x_{2}, 0 \leq x_{3}
\end{aligned}
$$

Show that the global minimum is $x^{*}=(1 / 2,1 / 2,0)$. Write a computer program implementing the conditional gradient method with the line minimization stepsize rule. (Here, there is a closed form expression for the minimizing stepsize.) Verify computationally that for a starting point $\left(\xi_{1}, \xi_{2}, \xi_{3}\right)$ with $\xi_{i}>0$ for all $i$, and $\xi_{1} \neq \xi_{2}$, the rate of convergence is not linear in the sense that

$$
\frac{f\left(x^{k+1}\right)-f\left(x^{*}\right)}{f\left(x^{k}\right)-f\left(x^{*}\right)} \rightarrow 1
$$

---

# 2.2 

Apply the conditional gradient method to the optimal routing problem of Example 1.3 in Section 2.1. Show that the feasible direction is generated by solving certain minimum length path problems.

## 2.3 (Sublinear Rate of Convergence of the Conditional Gradient Method)

Consider the two-dimensional problem

$$
\begin{aligned}
& \operatorname{minimize} f(x)=x_{1}^{2}+\left(x_{2}+1\right)^{2} \\
& \text { subject to }-1 \leq x_{1} \leq 1, \quad 0 \leq x_{2}
\end{aligned}
$$

with global minimum $x^{*}=(0,0)$. Successive iterates of the conditional gradient method are shown in Fig. 2.2.5. Show that

$$
\lim _{k \rightarrow \infty} \frac{f\left(x^{k+1}\right)-f\left(x^{*}\right)}{f\left(x^{k}\right)-f\left(x^{*}\right)}=1
$$

and therefore the rate of convergence of $\left\{f\left(x^{k}\right)-f\left(x^{*}\right)\right\}$ is not linear. Hint: This is a tricky problem. Show that for all $k \geq 1$, we have $f\left(x^{k}\right)-f\left(x^{*}\right)=$ $2\left(\cos \beta^{k}\right)^{2}-1=\sin \left(2 \alpha^{k}\right)$, where $\alpha^{k}$ and $\beta^{k}$ are the angles shown in Fig. 2.2.5.
![[ch_2_p27_img11.jpeg]]

Figure 2.2.5. Successive iterates of the conditional gradient method for the problem of Exercise 2.3. The directions used (asymptotically) become horizontal and the convergence rate is sublinear.

## 2.4 (Another Convergence Proof of the Conditional Gradient Method)

Assume that the gradient of $f$ satisfies

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|, \quad \forall x, y \in X
$$

---

where $L$ is a positive constant.
(a) Show that if $d$ is a descent direction at $x$, then

$$
\min _{\alpha \in[0,1]} f(x+\alpha d) \leq f(x)+\delta
$$

where $\delta$ is the negative scalar given by

$$
\delta= \begin{cases}\frac{1}{2} \nabla f(x)^{\prime} d, & \text { if } \nabla f(x)^{\prime} d+L\|d\|^{2}<0 \\ -\frac{\left|\nabla f(x)^{\prime} d\right|^{2}}{2 L R^{2}}, & \text { otherwise }\end{cases}
$$

where $R$ is the diameter of $X$, that is,

$$
R=\max _{x, y \in X}\|x-y\|
$$

Hint: Use the descent lemma (Prop. A. 24 of Appendix A), which shows that

$$
f(x+\alpha d) \leq f(x)+\alpha \nabla f(x)^{\prime} d+\frac{\alpha^{2} L}{2}\|d\|^{2}
$$

Minimize over $\alpha \in[0,1]$ both sides of this inequality.
(b) Consider the conditional gradient method

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}
$$

with the limited minimization stepsize rule. Show that every limit of the sequence $\left\{x^{k}\right\}$ is a stationary point. Hint: Argue that, if $\left\{x^{k}\right\}$ has a limit point and $\delta^{k}$ corresponds to $x^{k}$ as in part (a), then $\delta^{k} \rightarrow 0$, and, therefore, also $\nabla f\left(x^{k}\right)^{\prime} d^{k} \rightarrow 0$. Take the limit in the relation $\nabla f\left(x^{k}\right)^{\prime} d^{k} \leq \nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)$ for all $x \in X$.

# 2.5 

Consider the conditional gradient method, assuming that the gradient of $f$ satisfies

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|, \quad \forall x, y \in X
$$

where $L$ is a positive constant. Show that if the stepsize $\alpha^{k}$ is given by

$$
\alpha^{k}=\min \left\{1, \frac{\nabla f\left(x^{k}\right)^{\prime}\left(x^{k}-\bar{x}^{k}\right)}{L\left\|x^{k}-\bar{x}^{k}\right\|^{2}}\right\}
$$

then every limit point of $\left\{x^{k}\right\}$ is a stationary point. Hint: Use the line of analysis of Exercise 2.4.

---

# 2.6 (Termination Criterion) 

Show that if $f$ is convex, then in the conditional gradient method we have for all $k$,

$$
f\left(x^{k}\right) \geq \min _{x \in X} f(x) \geq f\left(x^{k}\right)+\nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right)
$$

and the upper and lower bound difference, $\nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right)$, tends to zero as $k \rightarrow \infty$.

### 2.7 (Simplicial Decomposition Method [Hol74], [Hoh77], [HLV87], [VeH93])

Assume that $f$ is convex and $X$ is a bounded polyhedron. Let $\bar{x}^{k}$ be defined as an extreme point of $X$ satisfying

$$
\bar{x}^{k}=\arg \min _{x \in X} \nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)
$$

and suppose that $x^{k+1}$ is a solution of the optimization problem

$$
\begin{aligned}
& \text { minimize } f(x) \\
& \text { subject to } x \in X^{k}
\end{aligned}
$$

where $X^{k}$ is the convex hull of $x^{0}$ and the extreme points $\bar{x}^{0}, \ldots, \bar{x}^{k}$,

$$
X^{k}=\left\{x \mid x=\beta x^{0}+\sum_{i=0}^{k} \gamma_{i} \bar{x}^{i}, \beta \geq 0, \gamma_{i} \geq 0, \beta+\sum_{i=0}^{k} \gamma_{i}=1\right\}
$$

Using the fact that the number of extreme points of $X$ is finite, show that the method finds a global minimum of $f$ over $X$ in a finite number of iterations.

### 2.8 (A Variation of the Simplicial Decomposition Method)

Consider the method of Exercise 2.7 with the difference that the set $X^{k}$ is any closed subset of the set

$$
\left\{x \mid x=\beta x_{0}+\sum_{i=0}^{k} \gamma_{i} \bar{x}^{i}, \beta \geq 0, \gamma_{i} \geq 0, \beta+\sum_{i=0}^{k} \gamma_{i}=1\right\}
$$

that contains the segment connecting $x^{k}$ and $\bar{x}^{k}$. Prove that every limit point of a sequence $\left\{x^{k}\right\}$ generated by this method is a stationary point. Hint: Use the fact that for every $k$, the cost reduction $f\left(x^{k}\right)-f\left(x^{k+1}\right)$ is larger than the cost reduction that would be obtained by using an Armijo rule along the segment connecting $x^{k}$ and $\bar{x}^{k}$.

---

# 2.9 (Zoutendijk's Method of Feasible Directions) 

Consider the problem

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } g_{j}(x) \leq 0, \quad j=1, \ldots, r
\end{aligned}
$$

where $f: \Re^{n} \mapsto \Re$ is continuously differentiable, and $g_{j}: \Re^{n} \mapsto \Re$ are convex and continuously differentiable functions. For any $\epsilon \geq 0$ and feasible vector $x$, define by $A(x ; \epsilon)$ the set of $\epsilon$-active constraints, that is,

$$
A(x ; \epsilon)=\left\{j \mid-\epsilon \leq g_{j}(x) \leq 0\right\}
$$

and let $\phi(x ; \epsilon)$ be the optimal value of the following linear programming problem in the vector $(d, z) \in \Re^{n+1}$

$$
\begin{array}{ll}
\operatorname{minimize} & z \\
\text { subject to } & \|d\|_{\infty} \leq 1, \quad \nabla f(x)^{\prime} d \leq z \\
& \nabla g_{j}(x)^{\prime} d \leq z, \quad \forall j \in A(x ; \epsilon)
\end{array}
$$

where $\|\cdot\|_{\infty}$ denotes the maximum norm.
Zoutendijk's method uses two scalars $\bar{\epsilon}>0$ and $\gamma \in(0,1)$, and determines a feasible descent direction $d^{k}$ at a nonstationary vector $x^{k}$ as follows:
$\left(d^{k}, \phi\left(x^{k} ; \epsilon^{k}\right)\right)$ is a solution of problem $\left(P_{x^{k}, \epsilon^{k}}\right)$, where $\epsilon^{k}=\gamma^{m_{k}} \bar{\epsilon}$ and $m_{k}$ is the first nonnegative integer $m$ for which

$$
\phi\left(x^{k}, \gamma^{m} \bar{\epsilon}\right) \leq-\gamma^{m} \bar{\epsilon}
$$

(a) Show that $d^{k}$ as defined above is obtained after solving a finite number of problems $\left(P_{x^{k}, \epsilon^{k}}\right)$. Furthermore, $d^{k}$ is a feasible descent direction at $x^{k}$.
(b) Prove that $\left\{d^{k}\right\}$ is gradient related, establishing stationarity of the limit points of $\left\{x^{k}\right\}$.

### 2.10 (Min-H Method for Optimal Control)

Consider the problem of finding sequences $u=\left(u_{0}, u_{1}, \ldots, u_{N-1}\right)$ and $x=$ $\left(x_{1}, x_{2}, \ldots, x_{N}\right)$ that minimize

$$
g_{N}\left(x_{N}\right)+\sum_{i=0}^{N-1} g_{i}\left(x_{i}, u_{i}\right)
$$

subject to the constraints

$$
u_{i} \in U_{i}, \quad i=0, \ldots, N-1
$$

---

$$
x_{i+1}=\phi_{i}\left(x_{i}\right)+B_{i} u_{i}, \quad i=0, \ldots, N-1, \quad x_{0}: \text { given }
$$

where $B_{i}$ are given matrices of appropriate dimension. We assume that the functions $g_{i}$ and $\phi_{i}$ are continuously differentiable, and the sets $U_{i}$ are convex. We further assume that $g_{i}\left(x_{i}, \cdot\right)$ is convex as a function of $u_{i}$ for all $x_{i}$. Consider the algorithm that, given $u^{k}=\left(u_{0}^{k}, \ldots, u_{N-1}^{k}\right)$, and the corresponding state and costate trajectories $x^{k}=\left(x_{1}^{k}, \ldots, x_{N}^{k}\right)$ and $p^{k}=\left(p_{1}^{k}, \ldots, p_{N}^{k}\right)$, respectively, generates $u^{k+1}$ by

$$
u_{i}^{k+1}=u^{k}+\alpha^{k}\left(\bar{u}^{k}-u_{i}^{k}\right), \quad i=0, \ldots, N-1
$$

where $\alpha^{k}$ is generated by the Armijo rule or the limited minimization rule, and $\bar{u}^{k}$ is given by

$$
\bar{u}_{i}^{k}=\arg \min _{u_{i} \in U_{i}} H_{i}\left(x_{i}^{k}, u_{i}, p_{i+1}^{k}\right), \quad i=0, \ldots, N-1
$$

with

$$
H_{i}\left(x_{i}, u_{i}, p_{i+1}\right)=g_{i}\left(x_{i}, u_{i}\right)+p_{i+1}^{\prime} f_{i}\left(x_{i}, u_{i}\right)
$$

denoting the Hamiltonian function. Show that this is a feasible direction method.

# 2.3 GRADIENT PROJECTION METHODS 

The conditional gradient method uses a feasible direction obtained by solving a subproblem with linear cost. Gradient projection methods use instead a subproblem with quadratic cost. While this subproblem may be more complex, the resulting convergence rate is typically better. There are many variations of gradient projection methods, which may be viewed collectively as the constrained analogs of the gradient methods of Section 2.2 .

### 2.3.1 Feasible Directions and Stepsize Rules Based on Projection

We first develop informally the main ideas regarding projection methods, leaving the detailed results and proofs for the next subsection.

The simplest gradient projection method is a feasible direction method of the form

$$
x^{k+1}=x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right)
$$

where

$$
\bar{x}^{k}=\left[x^{k}-s^{k} \nabla f\left(x^{k}\right)\right]^{+}
$$

Here, $[\cdot]^{+}$denotes projection on the set $X, \alpha^{k} \in(0,1]$ is a stepsize, and $s^{k}$ is a positive scalar. Thus, to obtain the vector $\bar{x}^{k}$, we take a step $-s^{k} \nabla f\left(x^{k}\right)$

---

along the negative gradient, as in [[Chapter 1#Steepest Descent|steepest descent]]. We then project the result $x^{k}-s^{k} \nabla f\left(x^{k}\right)$ on $X$, thereby maintaining feasibility.

We may also view the scalar $s^{k}$ as a stepsize. This becomes evident when we select $\alpha^{k}=1$ for all $k$, in which case $x^{k+1}=\bar{x}^{k}$ and the method becomes

$$
x^{k+1}=\left[x^{k}-s^{k} \nabla f\left(x^{k}\right)\right]^{+}
$$

(see Fig. 2.3.1). If $x^{k}-s^{k} \nabla f\left(x^{k}\right)$ is feasible, the gradient projection iteration becomes an unconstrained steepest descent iteration, in either the form (3.1) or the form (3.3). In this case, Eq. (3.1) takes the form $x^{k+1}=$ $x^{k}-\alpha^{k} \nabla f\left(x^{k}\right)$, while Eq. (3.3) takes the form $x^{k+1}=x^{k}-s^{k} \nabla f\left(x^{k}\right)$.

Note that we have $x^{*}=\left[x^{*}-s \nabla f\left(x^{*}\right)\right]^{+}$for all $s>0$ if and only if $x^{*}$ is stationary (see Fig. 2.3.2), so the method stops at such a point.

In order for the method to make practical sense, it is necessary that the projection operation is fairly simple. This will be so if $X$ has a relatively simple structure. For example, when the constraints are bounds on the variables,

$$
X=\left\{x \mid \alpha_{i} \leq x_{i} \leq \beta_{i}, i=1, \ldots, n\right\}
$$

the $i$ th coordinate of the projection of a vector $x$ is given by

$$
[x]_{i}^{+}= \begin{cases}\alpha_{i} & \text { if } x_{i} \leq \alpha_{i} \\ \beta_{i} & \text { if } x_{i} \geq \beta_{i} \\ x_{i} & \text { otherwise }\end{cases}
$$

![[ch_2_p32_img12.jpeg]]

Figure 2.3.1. Illustration of a few iterations of the gradient projection method for the case where $\alpha^{k}=1$ for all $k$. Note that when $x^{k}-s^{k} \nabla f\left(x^{k}\right)$ belongs to $X$, the iteration reduces to an unconstrained steepest descent iteration.

---

![[ch_2_p33_img13.jpeg]]

Figure 2.3.2. Illustrating why the gradient projection method stops if and only if it finds a stationary point $x^{*}$. By definition, $x^{*}$ is a stationary point if

$$
\nabla f\left(x^{*}\right)^{\prime}\left(x-x^{*}\right) \geq 0, \quad \forall x \in X
$$

which is equivalent to

$$
\left(\left(x^{*}-s \nabla f\left(x^{*}\right)\right)-x^{*}\right)^{\prime}\left(x-x^{*}\right) \leq 0, \quad \forall x \in X, s>0
$$

This holds if and only if $x^{*}$ is the projection of $x^{*}-s \nabla f\left(x^{*}\right)$ on $X$ (Prop. 2.1.3).

# Stepsize Selection and Convergence 

There are several stepsize selection procedures in the gradient projection method. The following two are straightforward, although we will later argue that they are often not the best choices.

## Limited Minimization Rule

Here

$$
s^{k}=s: \text { constant }, \quad k=0,1, \ldots
$$

and $\alpha^{k}$ is chosen by minimization over $[0,1]$, that is,

$$
f\left(x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right)\right)=\min _{\alpha \in[0,1]} f\left(x^{k}+\alpha\left(\bar{x}^{k}-x^{k}\right)\right)
$$

---

# Armijo Rule Along the Feasible Direction 

Here

$$
s^{k}=s: \text { constant }, \quad k=0,1, \ldots
$$

and $\alpha^{k}$ is chosen by the Armijo rule over the interval $[0,1]$. In particular, fixed scalars $\beta$ and $\sigma>0$, with $\beta \in(0,1)$ and $\sigma \in(0,1)$ are chosen, and we set $\alpha^{k}=\beta^{m_{k}}$, where $m_{k}$ is the first nonnegative integer $m$ for which

$$
f\left(x^{k}\right)-f\left(x^{k}+\beta^{m}\left(\bar{x}^{k}-x^{k}\right)\right) \geq-\sigma \beta^{m} \nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right)
$$

Our general convergence result for feasible direction methods (Prop. 2.2.1) applies to the gradient projection method with the above two stepsize rules, provided we can show that the direction sequence is gradient related. Indeed this can be proved (see the following Prop. 2.3.1), so every limit point of a sequence generated by the gradient projection method with the line minimization rule or the Armijo rule is stationary.

## Armijo Rule Along the Projection Arc

Here the stepsize $\alpha^{k}$ is fixed at unity,

$$
\alpha^{k}=1, \quad k=0,1, \ldots
$$

and the stepsize $s^{k}$ is determined by successive reduction until an Armijo-like inequality is satisfied. This means that $x^{k+1}$ is determined by an Armijo-like search on the "projection arc"

$$
\left\{x^{k}(s) \mid s>0\right\}
$$

where, for all $s>0, x^{k}(s)$ is defined by

$$
x^{k}(s)=\left[x^{k}-s \nabla f\left(x^{k}\right)\right]^{+}
$$

In particular, fixed scalars $\bar{s}, \beta$, and $\sigma>0$, with $\bar{s}>0, \beta \in(0,1)$, and $\sigma \in(0,1)$ are chosen, and we set $s^{k}=\beta^{m_{k}} \bar{s}$, where $m_{k}$ is the first nonnegative integer $m$ for which

$$
f\left(x^{k}\right)-f\left(x^{k}\left(\beta^{m} \bar{s}\right)\right) \geq \sigma \nabla f\left(x^{k}\right)^{\prime}\left(x^{k}-x^{k}\left(\beta^{m} \bar{s}\right)\right)
$$

Figure 2.3.3 illustrates the Armijo rule on a projection arc.
The convergence properties of the gradient projection method for the Armijo rules along the feasible direction and along the projection arc are essentially the same, although the convergence proofs differ substantially [see the following Prop. 2.3.3, which also shows that the stepsize rules are well defined, that is, a stepsize will be found after a finite number of trials based on the test (3.7)].

---

![[ch_2_p35_img14.jpeg]]

Figure 2.3.3. (a) Illustration of the successive points tested by the Armijo rule along the projection arc. (b) Illustration of the successive tests of the Armijo inequality (3.7).

# Constant Stepsize 

Here, $s^{k}$ is also fixed at some constant $s>0$ and $\alpha^{k}$ is fixed at unity, that is,

$$
s^{k}=s: \text { constant }, \quad \alpha^{k}=1, \quad k=0,1, \ldots
$$

Similar to unconstrained gradient methods, it is possible to show that the limit points of a sequence generated by the gradient projection method with a constant stepsize are stationary, provided $s$ is sufficiently small and the gradient satisfies a Lipschitz continuity condition (see the following Prop. 2.3.2).

## Diminishing Stepsize

Here, $\alpha^{k}$ is fixed at unity and

$$
s^{k} \rightarrow 0, \quad \sum_{k=0}^{\infty} s^{k}=\infty
$$

---

The convergence properties of the diminishing stepsize rule are similar to the ones of their unconstrained counterpart. We cannot guarantee descent at each iteration, but descent is more likely as the stepsize diminishes. The convergence rate tends to be slow, so this stepsize rule is used primarily for singular problems or when there are errors in the gradient calculation. Note that if the constraint set is bounded, the generated sequence $\{x^{k}\}$ is bounded, and this tends to enhance the convergence properties of the method (see Prop. 1.2.4 in Section 1.2). For this reason, an artificial bounded constraint set is sometimes introduced even for unconstrained problems, when a diminishing stepsize rule is used.

**Identification of the Set of Active Constraints***

The main difference between the Armijo rules along the feasible direction and along the projection arc, is that the iterates produced by the latter are more likely to be at the boundary of the constraint set $X$ than the iterates produced by the former. To formalize this idea, let us assume that $X$ is specified by linear inequality constraints

$$X = \{ x \mid a'_{j} x \leq b_{j}, j = 1, \ldots, r \} $$

The set of active constraints at a point $x \in X$ is the set of indices of the constraints that are satisfied as equations at $x$, that is, the set

$$A(x) = \{ j \mid a'_{j} x = b_{j}, j = 1, \ldots, r \} $$

We say that at the $k$th iteration the $j$th constraint is active, if $a'_{j} x^k = b_{j}$. It can be shown under fairly mild assumptions, that the set of active constraints at $x^*$ are finitely identified using the Armijo rule along the projection arc (this was shown in [Ber76c] for the case of bound constraints; for analysis of more general cases, see [GaB84], [Dun87], [BuM88], [Bon89a], [AlK90], [Fla92], [DeT93], [Wri93a]). By this we mean that the set of active constraints at $x^k$, is the same as the set of active constraints at $x^*$ for all sufficiently large $k$. By contrast an analogous result cannot be shown for the gradient projection method using the Armijo rule along the feasible direction. Generally, the property of finite identification of the active constraints is significant, because it leads to a sharper rate of convergence analysis and because it facilitates the fruitful combination of the gradient projection method with other iterations, which are applicable when the set of active constraints stays fixed.

**Rate of Convergence**

The convergence rate properties of the gradient projection method are essentially the same as those of the unconstrained steepest descent method. As an example, assume that $f$ is quadratic of the form

$$f(x) = \frac{1}{2} x'Qx - b'x,$$

---

where $Q$ is positive definite, and let $x^{*}$ denote the unique minimum of $f$ over $X$. Consider the case of a constant stepsize $\left(a^{k}=1\right.$ and $s^{k}=s$ for all $k$ ). Using the nonexpansive property of the projection [Prop. 2.1.3(c)] and the gradient formula $\nabla f(x)=Q x-b$, we have

$$
\begin{aligned}
\left\|x^{k+1}-x^{*}\right\| & =\left\|\left[x^{k}-s \nabla f\left(x^{k}\right)\right]^{+}-\left[x^{*}-s \nabla f\left(x^{*}\right)\right]^{+}\right\| \\
& \leq\left\|\left(x^{k}-s \nabla f\left(x^{k}\right)\right)-\left(x^{*}-s \nabla f\left(x^{*}\right)\right)\right\| \\
& =\left\|(I-s Q)\left(x^{k}-x^{*}\right)\right\| \\
& \leq \max \{|1-s m|,|1-s M|\}\left\|x^{k}-x^{*}\right\|
\end{aligned}
$$

where $m$ and $M$ are the minimum and maximum eigenvalues of $Q$. This is precisely the rate of convergence estimate obtained for the unconstrained steepest descent method with constant stepsize [cf. Eq. (3.1) and Fig. 1.3.1 in Section 1.3].

We conclude that the gradient projection method suffers from the same type of slow convergence as steepest descent. This motivates us to consider scaling, that is, application of the gradient projection method in a different coordinate system.

# Scaled Gradient Projection 

Scaling for gradient projection is derived similar to steepest descent. In particular, at the $k$ th iteration, let $H^{k}$ be a positive definite matrix and consider a transformation of variables defined by

$$
x=\left(H^{k}\right)^{-1 / 2} y
$$

Then, in the space of $y$, the problem is written as

$$
\begin{aligned}
& \text { minimize } h^{k}(y) \equiv f\left(\left(H^{k}\right)^{-1 / 2} y\right) \\
& \text { subject to } y \in Y^{k}
\end{aligned}
$$

where $Y^{k}$ is the set

$$
Y^{k}=\left\{y \mid\left(H^{k}\right)^{-1 / 2} y \in X\right\}
$$

The gradient projection iteration for this problem takes the form

$$
y^{k+1}=y^{k}+\alpha^{k}\left(\bar{y}^{k}-y^{k}\right)
$$

where

$$
\bar{y}^{k}=\left[y^{k}-s^{k} \nabla h^{k}\left(y^{k}\right)\right]^{+}
$$

Equivalently, $\bar{y}^{k}$ can be defined as the vector that minimizes the expression

$$
\left\|y-y^{k}+s^{k} \nabla h^{k}\left(y^{k}\right)\right\|^{2}=\left(s^{k}\right)^{2}\left\|\nabla h^{k}\left(y^{k}\right)\right\|^{2}+2 s^{k} \nabla h^{k}\left(y^{k}\right)^{\prime}\left(y-y^{k}\right)+\left\|y-y^{k}\right\|^{2}
$$

---

over $y \in Y^{k}$. By neglecting the constant term $\left(s^{k}\right)^{2}\left\|\nabla h^{k}\left(y^{k}\right)\right\|^{2}$ and by dividing by $2 s^{k}$ in the right-hand side of the above expression, we see that

$$
\bar{y}^{k}=\arg \min _{y \in Y^{k}}\left\{\nabla h^{k}\left(y^{k}\right)^{\prime}\left(y-y^{k}\right)+\frac{1}{2 s^{k}}\left\|y-y^{k}\right\|^{2}\right\}
$$

By using the identifications

$$
\begin{gathered}
x=\left(H^{k}\right)^{-1 / 2} y, \quad x^{k}=\left(H^{k}\right)^{-1 / 2} y^{k}, \quad \bar{x}^{k}=\left(H^{k}\right)^{-1 / 2} \bar{y}^{k} \\
\nabla h^{k}\left(y^{k}\right)=\left(H^{k}\right)^{-1 / 2} \nabla f\left(x^{k}\right)
\end{gathered}
$$

and the definition (3.9) of the set $Y^{k}$, we see that the iteration (3.10) can be written as

$$
x^{k+1}=x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right)
$$

where $\bar{x}^{k}$ is given by [cf. Eq. (3.11)]

$$
\bar{x}^{k}=\arg \min _{x \in X}\left\{\nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)+\frac{1}{2 s^{k}}\left(x-x^{k}\right)^{\prime} H^{k}\left(x-x^{k}\right)\right\}
$$

We refer to the iteration defined by the two preceding equations as the scaled gradient projection method.

Note that we can view the quadratic problem in Eq. (3.13) as a generalized projection problem. In particular, it can be verified that $\bar{x}^{k}$ is the vector of $X$, which is at minimum distance from the vector $x^{k}-$ $s^{k}\left(H^{k}\right)^{-1} \nabla f\left(x^{k}\right)$, but with distance measured in terms of the norm $\|z\|_{H^{k}}=$ $\sqrt{x^{\prime} H^{k} z}$ (see Exercise 3.1).

Note also that positive definiteness of the scaling matrix $H^{k}$ is not strictly necessary for the validity of the method. What is needed is that $H^{k}$ satisfies

$$
\left(x-x^{k}\right)^{\prime} H^{k}\left(x-x^{k}\right)>0, \quad \forall x \in X \text { with } x \neq x^{k}
$$

In particular, if the set $X$ is contained in some linear manifold of $\Re^{n}, H^{k}$ need only be positive definite over the subspace that is parallel to this manifold.

The convergence properties of the scaled gradient projection method are the same as the ones of the unscaled version, provided that the sequence $\left\{H^{k}\right\}$ satisfies a condition guaranteeing that the generated direction sequence is gradient related (see the following Prop. 2.3.4).

The convergence rate of the scaled gradient projection method is governed by $m^{k}$ and $M^{k}$, the smallest and largest eigenvalues of the Hessian $\nabla^{2} H^{k}\left(y^{k}\right)$, which is equal to $\left(H^{k}\right)^{-1 / 2} \nabla^{2} f\left(x^{k}\right)\left(H^{k}\right)^{-1 / 2}$. As in unconstrained optimization, this suggests that one should try to choose the scaling matrix $H^{k}$ as close as possible to the Hessian matrix $\nabla^{2} f\left(x^{k}\right)$. Using a diagonal approximation to the Hessian is a particularly useful choice because it maintains the simplicity of the quadratic subproblem of Eq. (3.13), when the constraint set $X$ is simple (see Exercise 3.3). If $H^{k}=\nabla^{2} f\left(x^{k}\right)$, we obtain a constrained version of Newton's method, which we proceed to describe.

---

# Constrained Newton's Method 

Let us assume that $f$ is twice continuously differentiable and that the Hessian matrix $\nabla^{2} f\left(x^{k}\right)$ is positive definite for all $x \in X$. Consider the scaled gradient projection method with scaling matrix $H^{k}=\nabla^{2} f\left(x^{k}\right)$. It is given by

$$
x^{k+1}=x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right)
$$

where

$$
\bar{x}^{k}=\arg \min _{x \in X}\left\{\nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)+\frac{1}{2 s^{k}}\left(x-x^{k}\right)^{\prime} \nabla^{2} f\left(x^{k}\right)\left(x-x^{k}\right)\right\}
$$

Note that if $s^{k}=1$, the quadratic cost above is the second order Taylor series expansion of $f$ around $x^{k}$ (except for a constant term). In particular, when

$$
\alpha^{k}=1, \quad s^{k}=1
$$

$x^{k+1}$ is the vector that minimizes the second order Taylor series expansion around $x^{k}$, just as in the case of unconstrained optimization. We expect, therefore, that for a starting point $x^{0}$ sufficiently close to a local minimum $x^{*}$, the method with unity stepsizes $\alpha^{k}$ and $s^{k}$ converges to $x^{*}$ superlinearly. Indeed this can be shown using a nearly identical proof to the one for the corresponding unconstrained case (see the following Prop. 2.3.5).

When a good starting point is unknown, one of the Armijo rules with unity initial stepsize may be used to improve the convergence properties of the method. The convergence and rate of convergence results that can be proved for such methods hold no surprises and are very similar to the corresponding unconstrained optimization results of Sections 1.3 and 1.4 (see Exercise 3.2).

The main difficulty with Newton's method is that the quadratic direction finding subproblem (3.15) may not be simple, even when the constraint set $X$ has a simple structure. Thus the method typically makes practical sense only for problems of small dimension. This motivates the development of methods that attain the fast convergence of Newton's method, while maintaining the simplicity of the direction finding subproblem when $X$ is a simple set. We will discuss methods of this type in Section 2.4.

## Variations of the Scaled Gradient Projection Method*

There are several variations of the scaled gradient projection method. We describe some of the simpler possibilities here and we discuss three major variations in Sections $2.4,2.5$, and 2.6. In all these variations, $X$ is specified by linear inequality constraints

$$
X=\left\{x \mid a_{j}^{\prime} x \leq b_{j}, j=1, \ldots, r\right\}
$$

---

# Projecting on an Expanded Constraint Set 

In this variation, we modify the direction finding quadratic subproblem of Eq. (3.13) by considering only the active or nearly active constraints at the current point $x^{k}$. In particular, we fix a positive scalar $\epsilon$ and obtain a vector $\tilde{x}^{k}$ as

$$
\tilde{x}^{k}=\arg \min _{x \in X_{k}}\left\{\nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)+\frac{1}{2 s^{k}}\left(x-x^{k}\right)^{\prime} H^{k}\left(x-x^{k}\right)\right\}
$$

where

$$
X^{k}=\left\{x \mid a_{j}^{\prime} x \leq b_{j}, \text { for all } j \text { with } b_{j}-\epsilon \leq a_{j}^{\prime} x^{k} \leq b_{j}\right\}
$$

Note that $X^{k}$ potentially involves a much smaller number of constraints than $X$, so there may be an advantage in solving the quadratic subproblem of Eq. (3.16) rather than the subproblem of Eq. (3.13). Note also that $\tilde{x}^{k}$ need not be a feasible point because $X^{k}$ may contain $X$ strictly. A feasible descent direction can be obtained, however, as $\bar{x}^{k}-x^{k}$, where

$$
\bar{x}^{k}=\bar{\gamma} \tilde{x}^{k}+(1-\bar{\gamma}) x^{k}
$$

and $\bar{\gamma}$ is the largest $\gamma \in[0,1]$ such that $\gamma \tilde{x}^{k}+(1-\gamma) x^{k} \in X$. The next iterate $x^{k+1}$ is obtained by the usual formula

$$
x^{k+1}=x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right)
$$

where $\alpha^{k}$ is calculated using either the limited minimization or the Armijo rule.

The convergence and rate of convergence analysis of this section apply with minor modifications to the preceding variation; see Exercise 3.4.

## Projecting on a Restricted Constraint Set

This variation is similar to the simplicial decomposition method described in Exercise 2.7 in Section 2.2.

At the $k$ th iteration, let $\tilde{x}^{k}$ be defined as an extreme point of $X$ satisfying

$$
\tilde{x}^{k}=\arg \min _{x \in X} \nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)
$$

and define $\bar{x}^{k}$ as

$$
\bar{x}^{k}=\arg \min _{x \in X^{k}}\left\{\nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)+\frac{1}{2 s^{k}}\left(x-x^{k}\right)^{\prime} H^{k}\left(x-x^{k}\right)\right\}
$$

where $X^{k}$ is the convex hull of $x^{0}$ and the extreme points $\tilde{x}^{0}, \ldots, \tilde{x}^{k}$,

$$
X^{k}=\left\{x \mid x=\beta x^{0}+\sum_{i=0}^{k} \gamma_{i} \tilde{x}^{i}, \beta \geq 0, \gamma_{i} \geq 0, \beta+\sum_{i=0}^{k} \gamma_{i}=1\right\}
$$

---

Note that $X^{k}$ is the convex hull of the preceding subset $X^{k-1}$ and the extreme point $\tilde{x}^{k}$ and that it is possible that $\tilde{x}^{k} \in X^{k-1}$, in which case $X^{k}=X^{k-1}$. Using the fact that the number of extreme points of $X$ is finite, we can show that the method

$$
x^{k+1}=x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right)
$$

with $\alpha^{k}$ obtained by either the limited minimization or the Armijo rule has the same convergence properties as the regular scaled gradient projection method (Exercise 3.5). Note that by using the subset $X^{k}$ in place of $X$, the direction finding subproblem may be greatly simplified. For an example where this simplification is important, see Exercise 3.6.

# Combinations with Unconstrained Optimization Methods 

We mentioned earlier that when the Armijo rule along the projection arc is used, the gradient projection method tends to identify the active constraints at a minimum in a finite number of iterations. Thus, for sufficiently large $k$ the method becomes equivalent to steepest descent over a linear manifold

$$
\left\{x \mid a_{j}^{\prime} x=b_{j}, j \in A\left(x^{*}\right)\right\}
$$

where $A\left(x^{*}\right)$ is the set of active constraints at the eventual limit $x^{*}$. Once this occurs, one may introduce iterations that implement a Newton-like method along this manifold. There is an endless range of possibilities here; any kind of unconstrained minimization method can be adapted fruitfully in this context. The stepsize rule for these methods should be modified so that the generated iterates are feasible, i.e. the stepsize $\alpha^{k}$ is small enough so that the vector $x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right)$ does not violate any of the inactive constraints $a_{j}^{\prime} x \leq b_{j}$, $j \notin A\left(x^{k}\right)$.

A difficulty with the type of combined method just described is that one can never really be sure that the active constraints have been identified. Thus a trial and error approach is usually followed, whereby one switches to the unconstrained method after one or more consecutive gradient projection iterations do not change the active set, and one switches to the gradient projection method following a change in the set of active constraints. There are many possibilities along these lines, some of which may be found in [Ber76c], [DeT83], [MoT89], [CGT91].

## The Main Limitation of the Gradient Projection Method - Alternatives

The principal drawback of the gradient projection method is the substantial overhead for computing the projection at each iteration. To obtain a good convergence rate, nondingmal (Newton-like) scaling must often be used, as for example in the iteration (3.11)-(3.15). Unfortunately, in this case, even if the constraints are simple, such as bounds on the variables, the corresponding quadratic program can be very time-consuming.

---

The next three sections discuss different ways to overcome the difficulty for the case where the constraint set is a polyhedron. This is done by solving equality constrained quadratic programs in place of inequality constrained quadratic programs. Solving an equality constrained quadratic program is much simpler, because it amounts to projection on a subspace and involves just the solution of an associated system of linear equations (see Example 1.5 in Section 2.1).

In Section 2.4, we discuss the two-metric projection method, given by

$$
x^{k+1}=\left[x^{k}-\alpha^{k} D^{k} \nabla f\left(x^{k}\right)\right]^{+}
$$

for the case of bound constraints. This is a natural and simple adaptation of unconstrained Newton-like methods. The main difficulty here is that an arbitrary positive definite matrix $D^{k}$ will not necessarily yield a descent direction. However, it turns out that if some of the off-diagonal terms of $D^{k}$ are zero, one can obtain descent. Furthermore, one can select $D^{k}$ as the inverse of a partially diagonalized version of the Hessian matrix $\nabla^{2} f\left(x^{k}\right)$ and attain the superlinear rate of Newton's method.

In Section 2.5, we discuss manifold suboptimization methods, which are based on scaled gradient projection on the manifold of active constraints. The motivation for these methods is similar to the one for the previously discussed combinations of gradient projection iterations (aimed at identifying the active constraints) with unconstrained iterations on the manifold of active constraints. The main difference is that at the typical iteration of a manifold suboptimization method, at most one constraint can be added or subtracted from the active set.

Finally, in Section 2.6, we discuss some interior point methods, which generate sequences of points that are interior to the constraint set. These methods have proved particularly effective for linear programming problems. They can be viewed as scaled projection methods, but the choice of the scaling matrix is crucial for good performance. We discuss the so called affine scaling method, which has performed well in practice.

# 2.3.2 Convergence Analysis* 

Generally, the convergence results for (unscaled) gradient projection are similar to those for steepest descent. The details, however, are somewhat more complicated, particularly for the Armijo rule along the projection arc. The following three propositions, each relating to different stepsize rules, are the main results.

Proposition 2.3.1: (Armijo Rule and Limited Minimization Rule Along the Descent Direction) Let $\left\{x^{k}\right\}$ be a sequence generated by the gradient projection method with $\alpha^{k}$ chosen by the limited

---

minimization rule or by the Armijo rule along the feasible direction. Then every limit point of $\left\{x^{k}\right\}$ is stationary.

Proof: We will show that the direction sequence $\left\{\bar{x}^{k}-x^{k}\right\}$ is gradient related. Then, application of Prop. 2.2.1 proves the result. Indeed, suppose that $\left\{x^{k}\right\}_{k \in K}$ converges to a nonstationary point $\bar{x}$. We must prove that

$$
\begin{gathered}
\limsup _{k \rightarrow \infty, k \in K}\left\|\bar{x}^{k}-x^{k}\right\|<\infty \\
\limsup _{k \rightarrow \infty, k \in K} \nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right)<0
\end{gathered}
$$

By the continuity of the projection [Prop. 2.1.3(c)], we have

$$
\lim _{k \rightarrow \infty, k \in K} \bar{x}^{k}=[\tilde{x}-s \nabla f(\tilde{x})]^{+}
$$

so Eq. (3.24) holds because $\left\{\left\|\bar{x}^{k}-x^{k}\right\|\right\}_{k \in K}$ converges to $\|[\tilde{x}-s \nabla f(\tilde{x})]^{+}-$ $\tilde{x} \|$. To show Eq. (3.25), we note that, by the characteristic property of the projection [Prop. 2.1.3(b)], we have

$$
\left(x^{k}-s \nabla f\left(x^{k}\right)-\bar{x}^{k}\right)^{\prime}\left(x-\bar{x}^{k}\right) \leq 0, \quad \forall x \in X
$$

Applying this relation with $x=x^{k}$, we obtain

$$
\nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right) \leq-\frac{1}{s}\left\|x^{k}-\bar{x}^{k}\right\|^{2}
$$

By taking limit in the above relation, we obtain

$$
\limsup _{k \rightarrow \infty, k \in K} \nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right) \leq-\frac{1}{s}\left\|\tilde{x}-[\tilde{x}-s \dot{\nabla} f(\tilde{x})]^{+}\right\|^{2}
$$

Since $\tilde{x}$ is nonstationary, the right-hand side of the above inequality is negative (cf. Fig. 2.3.2), proving Eq. (3.25). Q.E.D.

Proposition 2.3.2: (Constant Stepsize) Let $\left\{x^{k}\right\}$ be a sequence generated by the gradient projection method with $\alpha^{k}=1$ and $s^{k}=s$ for all $k$. Assume that for some constant $L>0$, we have

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|, \quad \forall x, y \in X
$$

Then, if $0<s<2 / L$, every limit point of $\left\{x^{k}\right\}$ is stationary.

---

Proof: By the using the descent lemma (Prop. A. 24 in Appendix A), we have

$$
f\left(x^{k+1}\right)-f\left(x^{k}\right)=f\left(\bar{x}^{k}\right)-f\left(x^{k}\right) \leq \nabla f\left(x^{k}\right)^{\prime}\left(\bar{x}^{k}-x^{k}\right)+\frac{L}{2}\left\|\bar{x}^{k}-x^{k}\right\|^{2}
$$

Using this equation and Eq. (3.27), we obtain

$$
f\left(x^{k+1}\right)-f\left(x^{k}\right) \leq\left(\frac{L}{2}-\frac{1}{s}\right)\left\|\bar{x}^{k}-x^{k}\right\|^{2}
$$

If $\left\{x^{k}\right\}$ has a limit point, the left-hand side of the above equation tends to zero, while if $s<2 / L$, the right-hand side is nonpositive. Therefore, $\left\|\bar{x}^{k}-x^{k}\right\| \rightarrow 0$, which implies that for every limit point $\tilde{x}$ of $\left\{x^{k}\right\}$ we have $[\tilde{x}-s \nabla f(\tilde{x})]^{+}=\tilde{x}$ [cf. Eq. (3.26)], so $\tilde{x}$ is stationary (cf. Fig. 2.3.2). Q.E.D.

An extension of the above proposition that substantially weakens the Lipschitz assumption (3.28), is given in Exercise 3.7. The next proposition, originally shown in [GaB82], requires a fairly sophisticated proof:

# Proposition 2.3.3: (Armijo Rule Along the Projection Arc) 

(a) For every $x \in X$ there exists a scalar $s_{x}>0$ such that

$$
f(x)-f(x(s)) \geq \sigma \nabla f(x)^{\prime}(x-x(s)), \quad \forall s \in\left[0, s_{x}\right]
$$

where $x(s)=[x-s \nabla f(x)]^{+}$[cf. Eq. (3.6)].
(b) Let $\left\{x^{k}\right\}$ be a sequence generated by the gradient projection method with $\alpha^{k}=1$ for all $k$ and with the stepsize $s^{k}$ chosen by the Armijo rule along the projection arc. Then every limit point of $\left\{x^{k}\right\}$ is stationary.

Proof: We first show the following lemma.

Lemma 2.3.1: For every $x \in X$ and $z \in \Re^{n}$, the function $g:[0, \infty) \mapsto$ $\Re$ defined by

$$
g(s)=\frac{\left\|[x+s z]^{+}-x\right\|}{s}, \quad \forall s>0
$$

is monotonically nonincreasing.

---

Proof: Fix $x \in X, z \in \Re^{n}$, and $\gamma>1$. Denote

$$
a=x+z, \quad b=x+\gamma z
$$

Let $\bar{a}$ and $\bar{b}$ be the projections on $X$ of $a$ and $b$, respectively. It will suffice to show that $g(\gamma) \leq g(1)$, or equivalently,

$$
\|\bar{b}-x\| \leq \gamma\|\bar{a}-x\|
$$

If $\bar{a}=x$ then clearly $\bar{b}=x$, so Eq. (3.33) holds. Also if $a \in X$, then $\bar{a}=a=x+z$, so Eq. (3.33) becomes $\|\bar{b}-x\| \leq \gamma\|z\|=\|b-x\|$, which again holds since the projection is a nonexpansive mapping [cf. Prop. 2.1.3(c)]. Finally, if $\bar{a}=\bar{b}$, then Eq. (3.33) also holds. Therefore, it will suffice to show Eq. (3.33) in the case where $\bar{a} \neq \bar{b}, \bar{a} \neq x, \bar{b} \neq x, a \notin X$, shown in Fig. 2.3.4.
![[ch_2_p45_img15.jpeg]]

Figure 2.3.4. Construction used in the proof of Lemma 2.3.1.

Let $H_{a}$ and $H_{b}$ be the two hyperplanes that are orthogonal to $\bar{b}-\bar{a}$ and pass through $\bar{a}$ and $\bar{b}$, respectively. Since $(\bar{b}-\bar{a})^{\prime}(b-\bar{b}) \geq 0$ and $(\bar{b}-\bar{a})^{\prime}(a-\bar{a}) \leq 0$, we have that neither $a$ nor $b$ lie strictly between the two hyperplanes $H_{a}$ and $H_{b}$. Furthermore, $x$ lies on the same side of $H_{a}$ as $a$, so $x \notin H_{a}$. Denote the intersections of the line $\{x+\alpha(b-x) \mid \alpha \in \Re\}$ with $H_{a}$ and $H_{b}$ by $s_{a}$ and $s_{b}$, respectively. Denote the intersection of the line $\{x+\alpha(\bar{a}-x) \mid \alpha \in \Re\}$ with $H_{b}$ by $w$. We have

$$
\begin{aligned}
\gamma & =\frac{\|b-x\|}{\|a-x\|} \geq \frac{\left\|s_{b}-x\right\|}{\left\|s_{a}-x\right\|} \\
& =\frac{\|w-x\|}{\|\bar{a}-x\|}=\frac{\|w-\bar{a}\|+\|\bar{a}-x\|}{\|\bar{a}-x\|} \\
& \geq \frac{\|\bar{b}-\bar{a}\|+\|\bar{a}-x\|}{\|\bar{a}-x\|} \geq \frac{\|\bar{b}-x\|}{\|\bar{a}-x\|}
\end{aligned}
$$

---

where the second equality is by similarity of triangles, the next to last inequality follows from the orthogonality relation $(w-\bar{b})^{\prime}(\bar{b}-\bar{a})=0$, and the last inequality is obtained from the triangle inequality. From Eq. (3.34), we obtain Eq. (3.33), which was to be proved. Q.E.D.

We now return to the proof of Prop. 2.3.3:
(a) By the projection theorem, we have

$$
(x-x(s))^{\prime}(x-s \nabla f(x)-x(s)) \leq 0, \quad \forall x \in X, s>0
$$

Hence

$$
\nabla f(x)^{\prime}(x-x(s)) \geq \frac{\|x-x(s)\|^{2}}{s} \quad \forall x \in X, s>0
$$

If $x$ is stationary, the conclusion holds with $s_{x}$ any positive scalar, so assume that $x$ is nonstationary and therefore $\|x-x(s)\| \neq 0$ for all $s>0$. By the mean value theorem, we have for all $x \in X$ and $s \geq 0$,

$$
f(x)-f(x(s))=\nabla f(x)^{\prime}(x-x(s))+\left(\nabla f\left(\xi_{s}\right)-\nabla f(x)\right)^{\prime}(x-x(s))
$$

where $\xi_{s}$ lies on the line segment joining $x$ and $x(s)$. Therefore, Eq. (3.30) can be written as

$$
(1-\sigma) \nabla f(x)^{\prime}(x-x(s)) \geq\left(\nabla f(x)-\nabla f\left(\xi_{s}\right)\right)^{\prime}(x-x(s))
$$

From Eq. (3.36) and Lemma 2.3.1, we have for all $s \in(0,1]$,

$$
\nabla f(x)^{\prime}(x-x(s)) \geq \frac{\|x-x(s)\|^{2}}{s} \geq\|x-x(1)\| \cdot\|x-x(s)\|
$$

Therefore, Eq. (3.37) is satisfied for all $s \in(0,1]$ such that

$$
(1-\sigma)\|x-x(1)\| \geq\left(\nabla f(x)-\nabla f\left(\xi_{s}\right)\right)^{\prime} \frac{x-x(s)}{\|x-x(s)\|}
$$

It is seen that there exists $s_{x}>0$ such that the above relation, and therefore also Eqs. (3.37) and (3.30), are satisfied for $s \in\left(0, s_{x}\right]$.
(b) Part (a), together with Eq. (3.36) and the definition of the Armijo rule along the projection arc [cf. Eqs. (3.6) and (3.7)], show that $s^{k}$ is well defined as a positive number for all $k$. Let $\bar{x}$ be a limit point of $\left\{x^{k}\right\}$ and let $\left\{x^{k}\right\}_{K}$ be a subsequence converging to $\bar{x}$. Since $\left\{f\left(x^{k}\right)\right\}$ is monotonically nonincreasing, we have $f\left(x^{k}\right) \rightarrow f(\bar{x})$. Consider two cases:

---

Case 1: $\liminf _{k \rightarrow \infty, k \in K} s^{k}>\hat{s}$ for some $\hat{s}>0$. Then, from Eq. (3.36) and Lemma 2.3.1, we have for all $k \in K$ that are sufficiently large

$$
\begin{aligned}
f\left(x^{k}\right)-f\left(x^{k+1}\right) & \geq \sigma \nabla f\left(x^{k}\right)^{\prime}\left(x^{k}-x^{k+1}\right) \\
& \geq \sigma \frac{\left\|x^{k}-x^{k+1}\right\|^{2}}{s^{k}} \\
& =\frac{\sigma s^{k}\left\|x^{k}-x^{k+1}\right\|^{2}}{\left(s^{k}\right)^{2}} \\
& \geq \frac{\sigma \hat{s}\left\|x^{k}-x^{k}(\bar{s})\right\|^{2}}{\bar{s}^{2}}
\end{aligned}
$$

where $\bar{s}$ is the initial stepsize of the Armijo rule. Taking limit as $k \rightarrow \infty$, $k \in K$, we obtain

$$
0 \geq \frac{\sigma \hat{s}\|\bar{x}-\bar{x}(\bar{s})\|^{2}}{\bar{s}^{2}}
$$

Hence $\bar{x}=\bar{x}(\bar{s})$ and it follows that $\bar{x}$ is stationary (cf. Fig. 2.3.2).
Case 2: $\liminf _{k \rightarrow \infty, k \in \bar{K}} s^{k}=0$. Then there exists a subsequence $\left\{s^{k}\right\}_{\bar{K}}$, $\bar{K} \subset K$, converging to zero. It follows that for all $k \in \bar{K}$, which are sufficiently large, the Armijo test (3.7) will be failed at least once (i.e., $m^{k} \geq 1$ ) and therefore

$$
f\left(x^{k}\right)-f\left(x^{k}\left(\beta^{-1} s^{k}\right)\right)<\sigma \nabla f\left(x^{k}\right)^{\prime}\left(x^{k}-x^{k}\left(\beta^{-1} s^{k}\right)\right)
$$

Furthermore, for all such $k \in \bar{K}, x^{k}$ cannot be stationary, since for stationary $x^{k}$, we have $s^{k}=\bar{s}$. Therefore,

$$
\left\|x^{k}-x^{k}\left(\beta^{-1} s^{k}\right)\right\|>0
$$

By the mean value theorem, we have

$$
\begin{aligned}
f\left(x^{k}\right)-f\left(x^{k}\left(\beta^{-1} s^{k}\right)\right)= & \nabla f\left(x^{k}\right)^{\prime}\left(x^{k}-x^{k}\left(\beta^{-1} s^{k}\right)\right) \\
& +\left(\nabla f\left(\xi^{k}\right)-\nabla f\left(x^{k}\right)\right)^{\prime}\left(x^{k}-x^{k}\left(\beta^{-1} s^{k}\right)\right)
\end{aligned}
$$

where $\xi^{k}$ lies in the line segment joining $x^{k}$ and $x^{k}\left(\beta^{-1} s^{k}\right)$. Combining Eqs. (3.38) and (3.40), we obtain for all $k \in \bar{K}$ that are sufficiently large

$$
(1-\sigma) \nabla f\left(x^{k}\right)^{\prime}\left(x^{k}-x^{k}\left(\beta^{-1} s^{k}\right)\right)<\left(\nabla f\left(x^{k}\right)-\nabla f\left(\xi^{k}\right)\right)^{\prime}\left(x^{k}-x^{k}\left(\beta^{-1} s^{k}\right)\right)
$$

Using Eq. (3.36) and Lemma 2.3.1, we obtain

$$
\begin{aligned}
\nabla f\left(x^{k}\right)^{\prime}\left(x^{k}-x^{k}\left(\beta^{-1} s^{k}\right)\right) & \geq \frac{\left\|x^{k}-x^{k}\left(\beta^{-1} s^{k}\right)\right\|^{2}}{\beta^{-1} s^{k}} \\
& \geq \frac{1}{\bar{s}}\left\|x^{k}-x^{k}(\bar{s})\right\| \cdot\left\|x^{k}-x^{k}\left(\beta^{-1} s^{k}\right)\right\|
\end{aligned}
$$

---

Combining the preceding two relations, and using the Schwartz inequality, we obtain for all $k \in \bar{K}$ that are sufficiently large

$$
\begin{aligned}
\frac{1-\sigma}{s}\left\|x^{k}-x^{k}(\bar{s})\right\| \cdot & \left\|x^{k}-x^{k}\left(\beta^{-1} s^{k}\right)\right\| \\
& <\left(\nabla f\left(x^{k}\right)-\nabla f\left(\xi^{k}\right)\right)^{\prime}\left(x^{k}-x^{k}\left(\beta^{-1} s^{k}\right)\right) \\
& \leq\left\|\nabla f\left(x^{k}\right)-\nabla f\left(\xi^{k}\right)\right\| \cdot\left\|x^{k}-x^{k}\left(\beta^{-1} s^{k}\right)\right\|
\end{aligned}
$$

Using Eqs. (3.39) and (3.43), we obtain

$$
\frac{1-\sigma}{\bar{s}}\left\|x^{k}-x^{k}(\bar{s})\right\|<\left\|\nabla f\left(x^{k}\right)-\nabla f\left(\xi^{k}\right)\right\|
$$

Since $s^{k} \rightarrow 0$ and $x^{k} \rightarrow \bar{x}$ as $k \rightarrow \infty, k \in \bar{K}$, it follows that $\xi^{k} \rightarrow \bar{x}$, as $k \rightarrow \infty, k \in \bar{K}$. Taking the limit in the above equation as $k \rightarrow \infty, k \in \bar{K}$, we obtain

$$
\|\bar{x}-\bar{x}(\bar{s})\| \leq 0
$$

Hence $\bar{x}=\bar{x}(\bar{s})$ and it follows that $\bar{x}$ is stationary. Q.E.D.

Proposition 2.3.4: (Convergence of Scaled Gradient Projection) Let $\left\{x^{k}\right\}$ be a sequence generated by the scaled gradient projection method with $\alpha^{k}$ chosen by the limited minimization rule or by the Armijo rule along the feasible direction. Assume that, for some positive scalars $c_{1}$ and $c_{2}$, the scaling matrices $H^{k}$ satisfy

$$
c_{1}\|z\|^{2} \leq z^{\prime} H^{k} z \leq c_{2}\|z\|^{2}, \quad \forall z \in \Re^{n}, k=0,1, \ldots
$$

Then every limit point of $\left\{x^{k}\right\}$ is stationary.

Proof: Almost identical to the proof of Prop. 2.3.1; left for the reader. Q.E.D.

Proposition 2.3.5: (Convergence of Newton's Method) Let $f$ be twice continuously differentiable with positive definite Hessian for all $x \in X$ and let $x^{*}$ be a local minimum of $f$ over $X$. There exists a $\delta>0$ such that if $\left\|x_{0}-x^{*}\right\|<\delta$, then the sequence $\left\{x^{k}\right\}$ generated by the constrained version of Newton's method [Eqs. (3.14) and (3.15)] with $\alpha^{k}=1$ and $s^{k}=1$ for all $k$, satisfies $\left\|x^{k}-x^{*}\right\|<\delta$ for all $k$ and $x^{k} \rightarrow x^{*}$. Furthermore, $\left\|x^{k}-x^{*}\right\|$ converges to zero superlinearly.

---

Proof: (Abbreviated) Denote $H^{k}=\nabla^{2} f\left(x^{k}\right)$ and consider the vector norms defined by $\|z\|_{H^{k}}=\sqrt{z^{\prime} H^{k} z}$ for $z \in \Re^{n}$, and the corresponding norms $\|A\|_{H^{k}}=\sup _{z \neq 0}\left(\|A z\|_{H^{k}} /\left.\|z\|_{H^{k}}\right)\right.$ for $n \times n$ matrices $A$. Let also $[z]^{+}$denote the projection of $z$ on $X$ with respect to the norm $\|\cdot\|_{H^{k}}$ (cf. Exercise 3.1). We have, using the result of Exercise 3.1 and the nonexpansive property of projections [Prop. 2.3.1(b) in Section 2.1, which can be generalized for the case of projections with respect to the norm $\|\cdot\|_{H^{k}}$ ],

$$
\begin{aligned}
\left\|x^{k+1}-x^{*}\right\|_{H^{k}} & =\left\|\left[x^{k}-\left(H^{k}\right)^{-1} \nabla f\left(x^{k}\right)\right]^{+}-x^{*}\right\|_{H^{k}} \\
& \leq\left\|x^{k}-\left(H^{k}\right)^{-1} \nabla f\left(x^{k}\right)-x^{*}\right\|_{H^{k}}
\end{aligned}
$$

Using this relation, we obtain as in the proof of Prop. 1.4.1 [cf. Eq. (4.13) of Section 1.4]

$$
\begin{aligned}
\left\|x^{k+1}-x^{*}\right\|_{H^{k}} \leq M & \left(\int_{0}^{1}\left\|\nabla^{2} f\left(x^{k}\right)-\nabla^{2} f\left(x^{*}+t\left(x^{k}-x^{*}\right)\right)\right\|_{H^{k}} d t\right) \\
& \cdot\left\|x^{k}-x^{*}\right\|_{H^{k}}
\end{aligned}
$$

for some $M>0$. By continuity of $\nabla^{2} f$, we can take $\delta$ sufficiently small to ensure that for $\left\|x^{k}-x^{*}\right\|<\delta$, the term under the integral sign is arbitrarily small. From this, the superlinear convergence of $\left\|x^{k}-x^{*}\right\|$ to zero follows. Q.E.D.

# E XERCISES 

### 3.1 (Scaled Gradient Projection Method)

Let $H^{k}$ be a positive definite matrix and define $\bar{x}^{k}$ by

$$
\bar{x}^{k}=\arg \min _{x \in X}\left\{\nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)+\frac{1}{2 s^{k}}\left(x-x^{k}\right)^{\prime} H^{k}\left(x-x^{k}\right)\right\}
$$

Show that $\bar{x}^{k}$ solves the problem

$$
\begin{aligned}
& \operatorname{minimize} \frac{1}{2}\left\|x-\left(x^{k}-s^{k}\left(H^{k}\right)^{-1} \nabla f\left(x^{k}\right)\right)\right\|_{H^{k}}^{2} \\
& \text { subject to } x \in X
\end{aligned}
$$

where $\|\cdot\|_{H^{k}}$ is the norm defined by $\|z\|_{H^{k}}=\sqrt{z^{\prime} H^{k} z}$.

---

# 3.2 

Let $f$ be twice continuously differentiable with positive definite Hessian for all $x \in X$. Consider a sequence $\left\{x^{k}\right\}$ generated by the constrained Newton method of Eqs. (3.14) and (3.15), where for each $k$, either $s^{k}=1$ and $\alpha^{k}$ is chosen by means of the Armijo rule along the feasible direction, or $\alpha^{k}=1$ and $s^{k}$ is chosen by the Armijo rule along the projection arc with initial stepsize $s=1$. Assume that $x^{k}$ converges to a local minimum $x^{*}$ with positive definite Hessian. Show that $\left\{\left\|x^{k}-x^{*}\right\|\right\}$ converges superlinearly and that there exists an integer $\bar{k} \geq 0$ such that $s^{k}=1$ for all $k \geq \bar{k}$.

### 3.3 (Gradient Projection with Diagonal Scaling)

Consider the gradient projection method with a diagonal scaling matrix having the terms $H_{i}^{k}$ along the diagonal.
(a) Show that for the case of bound constraints, $X=\left\{x \mid \alpha_{i} \leq x_{i} \leq \beta_{i}, i=\right.$ $1, \ldots, n\}$, the $i$ th coordinate of the vector $\bar{x}^{k}$ is given by

$$
\bar{x}_{i}^{k}= \begin{cases}\alpha_{i} & \text { if } x_{i}^{k}-\frac{s^{k}}{H_{i}^{k}} \frac{\partial f\left(x^{k}\right)}{\partial x_{i}} \leq \alpha_{i} \\ \beta_{i} & \text { if } x_{i}^{k}-\frac{s^{k}}{H_{i}^{k}} \frac{\partial f\left(x^{k}\right)}{\partial x_{i}} \geq \beta_{i} \\ x_{i}^{k}-\frac{s^{k}}{H_{i}^{k}} \frac{\partial f\left(x^{k}\right)}{\partial x_{i}} & \text { otherwise }\end{cases}
$$

(b) Derive the form of $\bar{x}^{k}$ for the case of a simplex constraint.

## 3.4

Consider the unscaled method with projection on an expanded constraint set [Eqs. (3.16)-(3.19)]. Show that:
(a) $\bar{x}^{k}-x^{k}$ is a feasible descent direction when $x^{k}$ is nonstationary.
(b) Every limit point of a sequence $\left\{x^{k}\right\}$ generated by the method is stationary.
(c) The convergence rate estimate of Eq. (3.8) holds for $x^{k}$ sufficiently close to $x^{*}$.

## 3.5

Consider the unscaled method with projection on a restricted constraint set [Eqs. (3.20)-(3.23)]. Show that the conclusions of parts (a)-(c) of Exercise 3.4 hold.

---

# 3.6 (Gradient Projection for Optimal Routing [BeT89], [BeG92], [LuT94b]) 

Describe the application of the diagonally scaled gradient projection method with projection on a restricted constraint set to the optimal routing problem of Example 1.3 in Section 2.1. Why is this method preferable to the ordinary diagonally scaled gradient projection method? Hint: The number of all possible paths may be astronomical but typically only a few of these paths carry positive flow.

## 3.7

Suppose that the Lipschitz condition

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|, \quad \forall x, y \in X
$$

[cf. Eq. (3.28)] is replaced by the following two conditions:
(i) For every bounded set $A \subset X$, there exists some constant $L$ such that

$$
\|\nabla f(x)-\nabla f(y)\| \leq L\|x-y\|, \quad \forall x, y \in A
$$

(ii) The set $\{x \in X \mid f(x) \leq c\}$ is bounded for every $c \in \Re$.

Show that the convergence result of Prop. 2.3.2 remains valid provided that the constant stepsize $s$ is allowed to depend on the choice of the initial vector $x_{0}$. Hint: Choose a stepsize that guarantees that $x^{k}$ stays within the level set $\left\{x \in X \mid f(x) \leq f\left(x_{0}\right)\right\}$.

## 3.8 (The Proximal Minimization Algorithm)

Consider the algorithm

$$
x^{k+1}=\arg \min _{x \in X}\left\{f(x)+\frac{1}{2 c^{k}}\left\|x-x^{k}\right\|^{2}\right\}
$$

for the case of the quadratic function $f(x)=(1 / 2) x^{\prime} Q x+b^{\prime} x$, and for a sequence of positive scalars $\left\{c^{k}\right\}$ such that $Q+\left(1 / c^{k}\right) I$ is positive definite for all $k$. (This algorithm will be reconsidered for the case where $f$ is a general convex function in Subsection 5.4.6.)
(a) Show that this algorithm is equivalent to the scaled gradient projection method (3.12)-(3.13) for the special choice of the scaling matrix $H^{k}=$ $Q+\left(1 / c^{k}\right) I$ and for stepsizes $s^{k}=1, \alpha^{k}=1$. Hint: Write $f$ as

$$
f(x)=\nabla f\left(x^{k}\right)^{\prime}\left(x-x^{k}\right)+\frac{1}{2}\left(x-x^{k}\right)^{\prime} Q\left(x-x^{k}\right)+\nabla f\left(x^{k}\right)^{\prime} x^{k}-\frac{1}{2} x^{k^{\prime}} Q x^{k}
$$

(b) Consider the following generalized version of the algorithm:

$$
x^{k+1}=x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right)
$$

---

where $\alpha^{k} \in(0,1]$,

$$
\bar{x}^{k}=\arg \min _{x \in X}\left\{f(x)+\frac{1}{2}\left(x-x^{k}\right)^{\prime} M^{k}\left(x-x^{k}\right)\right\}
$$

and $M^{k}$ is such that $Q+M^{k}$ is positive definite. Show that it can be viewed as a scaled gradient projection method.
(c) Assume that $X=\Re^{n}, M^{k}=Q$, and that $\alpha^{k}=1 / 2$. Show that the algorithm of part (b) converges in a single step.

# 2.4 TWO-METRIC PROJECTION METHODS 

We mentioned in the preceding section that gradient projection methods make practical sense only when the projection can be carried out fairly easily. A typical example is the case of an unscaled or diagonally scaled gradient projection method and an orthant constraint

$$
X=\{x \mid x \geq 0\}
$$

On the other hand we have seen that the convergence rate of the diagonally scaled gradient projection method is often unacceptably slow. This motivated the use of more general, nondiagonal scaling. Then, however, the scaled projection becomes a fairly complex quadratic programming problem, even for the case of the simple orthant constraint (4.1).

To resolve the slow convergence/complex implementation dilemma one may think of the method

$$
x^{k+1}=\left[x^{k}-\alpha^{k} D^{k} \nabla f\left(x^{k}\right)\right]^{+}
$$

Here, $D^{k}$ is a positive definite, not necessarily diagonal, matrix, and $[\cdot]^{+}$ denotes the usual easy projection on the orthant with respect to the Euclidean norm. Thus, except for the projection operation, this method is identical with the usual scaled gradient iteration for unconstrained minimization. We call this method a two-metric projection method because, in contrast with the methods of the preceding section, it embodies two different scaling matrices: one is the matrix $D^{k}$ scaling the gradient and the other is the identity matrix, which is used in the (Euclidean) projection norm. This allows the possibility of using Hessian information to select $D^{k}\left[\text { e.g. } D^{k}=\left(\nabla^{2} f\left(x^{k}\right)\right)^{-1}\right.$, as in unconstrained Newton's method], while maintaining the simplicity of the Euclidean norm projection on the orthant.

---

In this section we discuss two-metric projection methods for the case of the orthant constraint (4.1), that is, the problem

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } x \geq 0
\end{aligned}
$$

A similar analysis applies to the case of a box constraint, that is, upper and lower bounds on the coordinates of $x$ (see Exercise 4.1), the case of a Cartesian product of simplices (see Exercise 4.2), and more general cases (see the references cited at the end of the chapter).

There is a fundamental difficulty with a two-metric gradient projection method: it is not in general a descent iteration. In particular we may have $f\left(x^{k+1}\right)>f\left(x^{k}\right)$ for all positive stepsizes $\alpha^{k}$. This situation is illustrated in Fig. 2.4.1, where it can be seen that with an unfavorable choice of $D^{k}$, the method does not even recognize (i.e., stop at) a stationary point.
![[ch_2_p53_img16.jpeg]]

Figure 2.4.1. (a) An example where the iteration

$$
x^{k+1}=\left[x^{k}-\alpha^{k} D^{k} \nabla f\left(x^{k}\right)\right]^{+}
$$

fails to make progress at a nonstationary point with a poor choice of the scaling matrix $D^{k}$. (b) An example where the iteration fails to stop at a stationary point with a poor choice of the scaling matrix $D^{k}$.

It turns out, however, that there is a class of nondiagonal matrices $D^{k}$ for which descent is guaranteed. This class is sufficiently wide to allow superlinear convergence when $D^{k}$ properly embodies second derivative information. We introduce this class of matrices and we prove its key properties.

---

Let us denote for all $x \geq 0$

$$
I^{+}(x)=\left\{i \mid x_{i}=0, \frac{\partial f(x)}{\partial x_{i}}>0\right\}
$$

We say that a symmetric $n \times n$ matrix $D$ with elements $d_{i j}$ is diagonal with respect to a subset of indices $I \subset\{1, \ldots, n\}$, if

$$
d_{i j}=0, \quad \forall i \in I, j=1, \ldots, n, j \neq i
$$

Proposition 2.4.1: Let $x \geq 0$ and let $D$ be a positive definite symmetric matrix which is diagonal with respect to $I^{+}(x)$. Denote

$$
x(\alpha)=[x-\alpha D \nabla f(x)]^{+}, \quad \forall \alpha \geq 0
$$

(a) The vector $x$ is a stationary point if and only if

$$
x=x(\alpha), \quad \forall \alpha \geq 0
$$

(b) If $x$ is not a stationary point, there exists a scalar $\bar{\alpha}>0$ such that

$$
f(x(\alpha))<f(x), \quad \forall \alpha \in(0, \bar{\alpha}]
$$

Proof: (a) By relabeling the coordinates of $x$ if necessary, we will assume that for some integer $r$ we have

$$
I^{+}(x)=\{r+1, \ldots, n\}
$$

The following proof also applies with minor modifications when $I^{+}(x)$ is empty. Since $D$ is diagonal with respect to $I^{+}(x)$, it has the form

$$
D=\left(\begin{array}{ccccc}
\bar{D} & 0 & 0 & \cdots & 0 \\
0 & d_{r+1} & 0 & \cdots & 0 \\
0 & 0 & d_{r+2} & \cdots & 0 \\
\vdots & \vdots & \vdots & \vdots & \vdots \\
0 & 0 & 0 & \cdots & d_{n}
\end{array}\right)
$$

where $\bar{D}$ is positive definite and $d_{i}>0, i=r+1, \ldots, n$. Denote

$$
p=D \nabla f(x)
$$

---

Assume that $x$ is a stationary point. Then, from the necessary optimality conditions [Eqs. (1.3) and (1.4) in Section 2.1] and the definition (4.3) of $I^{+}(x)$, we have

$$
\begin{gathered}
\frac{\partial f(x)}{\partial x_{i}}=0, \quad \forall i=1, \ldots, r \\
\frac{\partial f(x)}{\partial x_{i}}>0, \quad \forall i=r+1, \ldots, n
\end{gathered}
$$

These relations and the positivity of $d_{i}$ imply that

$$
\begin{gathered}
p_{i}=0, \quad \forall i=1, \ldots, r \\
p_{i}>0, \quad \forall i=r+1, \ldots, n
\end{gathered}
$$

Since $x_{i}(\alpha)=\left[x_{i}-\alpha p_{i}\right]^{+}$and $x_{i}=0$ for $i=r+1, \ldots, n$, it follows that $x_{i}(\alpha)=x_{i}$ for all $i$, and $\alpha \geq 0$.

Conversely assume that $x=x(\alpha)$ for all $\alpha \geq 0$. Then, we must have

$$
\begin{array}{ll}
p_{i}=0, & \forall i=1, \ldots, n \text { with } x_{i}>0 \\
p_{i} \geq 0, & \forall i=1, \ldots, n \text { with } x_{i}=0
\end{array}
$$

Now by the definition of $I^{+}(x)$, we have that if $x_{i}=0$ and $i \notin I^{+}(x)$, then $\partial f(x) / \partial x_{i} \leq 0$. This, together with the relations above, imply that

$$
\sum_{i=1}^{r} p_{i} \frac{\partial f(x)}{\partial x_{i}} \leq 0
$$

Since by Eqs. (4.7) and (4.8),

$$
\left[\begin{array}{c}
p_{1} \\
\vdots \\
p_{r}
\end{array}\right]=\bar{D}\left[\begin{array}{c}
\frac{\partial f(x)}{\partial x_{1}} \\
\vdots \\
\frac{\partial f(x)}{\partial x_{r}}
\end{array}\right]
$$

while $\bar{D}$ is positive definite, we also have

$$
\sum_{i=1}^{r} p_{i} \frac{\partial f(x)}{\partial x_{i}} \geq 0
$$

and it follows that

$$
p_{i}=\frac{\partial f(x)}{\partial x_{i}}=0, \quad \forall i=1, \ldots, r
$$

Since for $i=r+1, \ldots, n$, we have $\partial f(x) / \partial x_{i}>0$ and $x_{i}=0$, we see that $x$ is a stationary point.

---

(b) Let $r$ be as in the proof of part (a). For $i=r+1, \ldots, n$, we have $\partial f(x) / \partial x_{i}>0, x_{i}=0$, and from Eqs. (4.7) and (4.8), we see that

$$
x_{i}=x_{i}(\alpha)=0, \quad \forall \alpha \geq 0, \quad i=r+1, \ldots, n
$$

Consider the set of indices

$$
\begin{gathered}
I_{1}=\left\{i \mid x_{i}>0 \text { or }\left(x_{i}=0 \text { and } p_{i}<0\right), i=\dot{1}, \ldots, r\right\} \\
I_{2}=\left\{i \mid x_{i}=0 \text { and } p_{i} \geq 0, i=1, \ldots, r\right\}
\end{gathered}
$$

Let

$$
\alpha_{1}=\sup \left\{\alpha \mid x_{i}-\alpha p_{i} \geq 0, \forall i \in I_{1}\right\}
$$

Note that, in view of the definition of $I_{1}, \alpha_{1}$ is either positive or $+\infty$. Define the vector $\bar{p}$ with coordinates

$$
\bar{p}_{i}= \begin{cases}p_{i} & \text { if } i \in I_{1} \\ 0 & \text { if } i \in I_{2} \text { or } i \in I^{+}(x)\end{cases}
$$

Using Eqs. (4.9)-(4.13), we have

$$
x(\alpha)=x-\alpha \bar{p}, \quad \forall \alpha \in\left(0 ; \alpha_{1}\right)
$$

In view of the definition of $I_{2}$ and $I^{+}(x)$, we obtain

$$
\frac{\partial f(x)}{\partial x_{i}} \leq 0, \quad \forall i \in I_{2}
$$

and hence

$$
\sum_{i \in I_{2}} \frac{\partial f(x)}{\partial x_{i}} p_{i} \leq 0
$$

Now using Eqs. (4.13) and (4.16), we have

$$
\nabla f(x)^{\prime} \bar{p}=\sum_{i \in I_{1}} \frac{\partial f(x)}{\partial x_{i}} p_{i} \geq \sum_{i=1}^{r} \frac{\partial f(x)}{\partial x_{i}} p_{i}
$$

Since $x$ is not a stationary point, by part (a) and Eq. (4.14), we must have $x \neq x(\alpha)$ for some $\alpha>0$ and hence also, in view of Eq. (4.9), $p_{i} \neq 0$ for some $i \in\{1, \ldots, r\}$. In view of the positive definiteness of $\bar{D}$, and Eqs. $(4.7)$ and $(4.8)$,

$$
\sum_{i=1}^{r} \frac{\partial f(x)}{\partial x_{i}} p_{i}>0
$$

It follows from Eq. (4.17) that

$$
\nabla f(x)^{\prime} \bar{p}>0
$$

---

Combining this relation with Eq. (4.14) and the fact $\alpha_{1}>0$, yields that $\bar{p}$ is a feasible descent direction at $x$ and there exists a scalar $\bar{\alpha}>0$ for which the desired relation (4.6) is satisfied. Q.E.D.

Based on Prop. 2.4.1 we conclude that to guarantee descent, the matrix $D^{k}$ in the iteration

$$
x^{k+1}=\left[x^{k}-\alpha^{k} D^{k} \nabla f\left(x^{k}\right)\right]^{+}
$$

should be chosen diagonal with respect to a subset of indices that contains

$$
I^{+}\left(x^{k}\right)=\left\{i \mid x_{i}^{k}=0, \frac{\partial f\left(x^{k}\right)}{\partial x_{i}}>0\right\}
$$

However, it turns out that to guarantee convergence one should implement the iteration more carefully. The reason is that the set $I^{+}\left(x^{k}\right)$ exhibits an undesirable discontinuity at the boundary of the constraint set, whereby given a sequence $\left\{x^{k}\right\}$ of interior points that converges to a boundary point $\bar{x}$ the set $I^{+}\left(x^{k}\right)$ may be strictly smaller than the set $I^{+}(\bar{x})$. This causes difficulties in proving convergence of the algorithm and may have an adverse effect on its rate of convergence. To bypass these difficulties one may add to the set $I^{+}\left(x^{k}\right)$ the indices of those variables $x_{i}^{k}$ that are "near" zero and satisfy $\partial f\left(x^{k}\right) / \partial x_{i}>0$. With such a modification and with a variation of the Armijo rule on the projection arc given in Section 2.3, one can prove a satisfactory convergence result. It is also possible to construct Newtonlike algorithms, where the nondiagonal portion of the matrix $D^{k}$ consists of the inverse of the Hessian submatrix corresponding to the indices not in $I^{+}\left(x^{k}\right)$, and to show a quadratic rate of convergence result under the appropriate assumptions (see the references cited at the end of the chapter for details).

# E XER C ISES 

## 4.1

Consider the problem

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } \alpha \leq x \leq \beta
\end{aligned}
$$

where $\alpha$ and $\beta$ are given vectors. Formulate and prove a result that parallels Prop. 2.4.1.

---

# 4.2 (Simplex Constraints [Ber82b]) 

Consider the problem

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } x \geq 0, \quad \sum_{i=1}^{n} x_{i}=1
\end{aligned}
$$

For a feasible vector $x$, let $\bar{i}$ be an index such that $x_{\bar{i}}>0$, and to simplify notation, assume without loss of generality that $\bar{i}=n$. Let $y=\left(x_{1}, \ldots, x_{n-1}\right)$, let

$$
h(y)=f\left(x_{1}, \ldots, x_{n-1}, 1-\sum_{i=1}^{n-1} x_{i}\right)
$$

and let

$$
I^{+}(y)=\left\{i \mid y_{i}=0, \frac{\partial h(y)}{\partial y_{i}}>0\right\}
$$

Suppose that $D$ is an $(n-1) \times(n-1)$ positive definite matrix, which is diagonal with respect to $I^{+}(y)$, and define

$$
\begin{gathered}
y(\alpha)=\left[y-\alpha D \nabla h(y)\right]^{+}, \quad \alpha \geq 0 \\
x(\alpha)=\left(1-\sum_{i=1}^{y(\alpha)} y_{i}(\alpha)\right)
\end{gathered}
$$

Show that:
(a) $x$ is a stationary point if and only if $x=x(\alpha)$ for all $\alpha \geq 0$.
(b) If $x$ is not stationary, there exists a scalar $\bar{\alpha}>0$ such that

$$
f(x(\alpha))<f(x), \quad \forall \alpha \in(0, \bar{\alpha}]
$$

(c) State a descent method of the Newton type in terms of the first and second derivatives of $f$.

### 2.5 MANIFOLD SUBOPTIMIZATION - QUADRATIC PROGRAMMING

The methods of this section are feasible direction methods for the linearly constrained problem

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } a_{j}^{\prime} x \leq b_{j}, \quad j=1, \ldots, r
\end{aligned}
$$

---

They may be viewed as variants of gradient projection methods, the main difference being that, to obtain a feasible descent direction, the gradient is projected on a linear manifold of active constraints rather than on the entire constraint set (see Fig. 2.5.1). This greatly simplifies the projection and provides the main advantage of the methods of this section. Once the set of active constraints at a solution is identified, the methods behave identically with (scaled) steepest descent methods. The early portion of the computation may be viewed as a systematic effort to identify the set of active constraints at a solution by searching through a sequence of successive manifolds typically differing by a single constraint. Thus the number of iterations required to identify the set of active constraints is typically at least as large as the number of constraints whose active/inactive status is different at the starting point than at the solution. It follows that the methods of this section are well suited only for problems with a relatively small number of constraints.
![[ch_2_p59_img17.jpeg]]

Figure 2.5.1. Illustration of the manifold suboptimization method. The method searches through a sequence of manifolds. Two successive manifolds typically differ by at most one constraint. In (a), the optimal solution is at a vertex, as in the case of linear programs. In (b), the optimal solution is in the interior of the constraint set.

Throughout this section, we we assume that at every feasible point $x$, the set of vectors

$$
\left\{a_{j} \mid j \in A(x)\right\}
$$

is linearly independent, where $A(x)$ is the set of indices of active constraints at $x$

$$
A(x)=\left\{j \mid a_{j}^{\prime} x=b_{j}, j=1, \ldots, r\right\}
$$

This assumption can be relaxed at the expense of some technical complications [basically enough indices have to be dropped from $A(x)$ so that the remaining vectors $a_{j}$ are linearly independent, but the indices dropped must be chosen carefully so that they define redundant constraints].

---

The iteration at a vector $x^{k}$ proceeds roughly as follows:
We try to find a feasible descent direction from the "enlarged" subspace

$$
S\left(x^{k}\right)=\left\{d \mid a_{j}^{\prime} d=0, j \in A(x)\right\}
$$

This is the subspace which is parallel to the manifold of active constraints, so a small movement along any $d \in S\left(x^{k}\right)$ does not change the set of active constraints.

There are two possibilities:
(a) A feasible descent direction $d^{k} \in S\left(x^{k}\right)$ is found, in which case the next point $x^{k+1}$ is given by

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}
$$

where $\alpha^{k}$ is obtained by some stepsize rule within the interval of stepsizes

$$
\left\{\alpha>0 \mid x^{k}+\alpha^{k} d^{k} \text { is feasible }\right\}
$$

(b) No feasible descent direction $d \in S\left(x^{k}\right)$ can be found because $x^{k}$ is stationary over the manifold $x^{k}+S\left(x^{k}\right)$ of active constraints.
In case (b), there are two possibilities: either $x^{k}$ is stationary over the entire constraint set

$$
\left\{x \mid a_{j}^{\prime} x \leq b_{j}, j=1, \ldots, r\right\}
$$

in which case the algorithm stops, or else one of the active constraints, say $\bar{j}$, is relaxed and a feasible descent direction belonging to the subspace

$$
\bar{S}\left(x^{k}\right)=\left\{d \mid a_{j}^{\prime} d=0, j \in A\left(x^{k}\right), j \neq \bar{j}\right\}
$$

is obtained (see Fig. 2.5.2). We will describe shortly how $\bar{j}$ is selected and how the corresponding direction is computed.

Central to manifold optimization methods are quadratic programming problems of the form

$$
\begin{aligned}
& \operatorname{minimize} \nabla f\left(x^{k}\right)^{\prime} d+\frac{1}{2} d^{\prime} H^{k} d \\
& \text { subject to } d \in S\left(x^{k}\right)=\left\{d \mid a_{j}^{\prime} d=0, j \in A\left(x^{k}\right)\right\}
\end{aligned}
$$

where $H^{k}$ is a symmetric positive definite matrix. Note that the solution of this problem can be viewed as a scaled projection of the gradient $\nabla f\left(x^{k}\right)$ on the subspace $S\left(x^{k}\right)$. The unique optimal solution is (see Example 1.5 in Section 2.1)

$$
\begin{gathered}
d^{k}=-\left(H^{k}\right)^{-1}\left(\nabla f\left(x^{k}\right)+A^{k^{\prime}} \mu\right) \\
\mu=-\left(A^{k}\left(H^{k}\right)^{-1} A^{k^{\prime}}\right)^{-1} A^{k}\left(H^{k}\right)^{-1} \nabla f\left(x^{k}\right)
\end{gathered}
$$

---

![[ch_2_p61_img18.jpeg]]

Figure 2.5.2. Illustration of the manifold suboptimization method when $x^{k}$ is stationary over the current manifold of active constraints. In (a), $x^{k}$ is stationary over the entire constraint set and the algorithm stops. In (b), one of the active constraints is relaxed and the search continues over a larger manifold.
where $A^{k}$ is the matrix that has as rows the vectors $a_{j}, j \in A\left(x^{k}\right)$. [It is actually sufficient that $H^{k}$ be positive definite on just the subspace $S\left(x^{k}\right)$; then there will again be a unique solution $d^{k}$, but Eqs. (5.9a) and (5.9b) must be appropriately modified.] Let us describe the typical iteration of the method:

We first note that since $d^{k}$ is the optimal solution of the quadratic programming problem (5.8), and since the vector $d=0$ is feasible for this problem, we must have

$$
\nabla f\left(x^{k}\right)^{\prime} d^{k}+\frac{1}{2} d^{\prime k} H^{k} d^{k} \leq 0
$$

If $d^{k} \neq 0$, we see that

$$
\nabla f\left(x^{k}\right)^{\prime} d^{k} \leq-\frac{1}{2} d^{\prime k} H^{k} d^{k}<0
$$

and $d^{k}$ is a feasible descent direction at $x^{k}$. The iteration is then given by

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}
$$

where $\alpha^{k}$ is obtained by some stepsize rule (e.g., Armijo, limited minimization, etc.) over the interval of stepsizes $\alpha$ for which $x^{k}+\alpha d^{k}$ is feasible, that is, the set

$$
\left\{\alpha>0 \mid a_{j}^{\prime}\left(x^{k}+\alpha d^{k}\right) \leq b_{j}, j \notin A\left(x^{k}\right)\right\}
$$

[Note that, by construction, we have $a_{j}^{\prime} d^{k}=0$ for all $j \in A\left(x^{k}\right)$, so by moving along the direction $d^{k}$, none of the active constraints will be violated or become inactive.]

---

On the other hand, if $d^{k}=0$, it follows from Eq. (5.9a) that $\nabla f\left(x^{k}\right)+$ $A^{k^{\prime}} \mu=0$, or equivalently, that for some scalars $\mu_{j}, j \in A\left(x^{k}\right)$, we have

$$
\nabla f\left(x^{k}\right)+\sum_{j \in A\left(x^{k}\right)} \mu_{j} a_{j}=0
$$

Now note that the feasible directions at $x^{k}$ are the vectors $d$ with $a_{j}^{\prime} d \leq 0$, for all $j \in A\left(x^{k}\right)$, so if $\mu_{j} \geq 0$ for all $j \in A\left(x^{k}\right)$, then we have

$$
\nabla f\left(x^{k}\right)^{\prime} d=-\sum_{j \in A\left(x^{k}\right)} \mu_{j} a_{j}^{\prime} d \geq 0
$$

for all feasible directions $d$, so that $x^{k}$ is stationary. Hence if $d^{k}=0$ and $x^{k}$ is nonstationary, we must have $\mu_{\bar{j}}<0$ for some index $\bar{j} \in A\left(x^{k}\right)$ (see Fig. 2.5.3 for an interpretation). Let $\bar{d}^{k}$ be the unique solution of the problem

$$
\begin{aligned}
& \text { minimize } \nabla f\left(x^{k}\right)^{\prime} d+\frac{1}{2} d^{\prime} \bar{H}^{k} d \\
& \text { subject to } d \in \bar{S}\left(x^{k}\right)=\left\{d \mid a_{j}^{\prime} d=0, j \in A\left(x^{k}\right), j \neq \bar{j}\right\}
\end{aligned}
$$

and $\bar{H}^{k}$ is a positive definite symmetric matrix. We claim that $\bar{d}^{k}$ is a feasible descent direction, so that the iteration takes the form

$$
x^{k+1}=x^{k}+\alpha^{k} \bar{d}^{k}
$$

where $\alpha^{k}$ is a stepsize chosen as earlier to maintain feasibility of $x^{k+1}$.
![[ch_2_p62_img19.jpeg]]

Figure 2.5.3. Dropping a constraint when $x^{k}$ is a stationary point on the manifold of active constraints. In this example, $x^{k}$ is a vertex and coincides with the manifold of active constraints. Therefore, $x^{k}$ is by default stationary on the manifold of active constraints, and one of the two constraints must be dropped. This is the constraint $a_{1}^{\prime} x=b_{1}$ for which the corresponding scalar $\mu_{1}$ is negative (notice that dropping the other constraint does not allow further descent).

---

To verify this, we first show that $\bar{d}^{k} \neq 0$. Indeed, if $\bar{d}^{k}=0$, then

$$
\nabla f\left(x^{k}\right)+\sum_{j \in A\left(x^{k}\right), j \neq \bar{j}} \bar{\mu}_{j} a_{j}=0
$$

for some $\bar{\mu}_{j}, j \in A\left(x^{k}\right), j \neq \bar{j}$. By subtracting Eqs. (5.12) and (5.15), we see that

$$
\mu_{\bar{j}} a_{\bar{j}}=\sum_{j \in A\left(x^{k}\right), j \neq \bar{j}}\left(\bar{\mu}_{j}-\mu_{j}\right) a_{j}
$$

which, in view of the choice $\mu_{\bar{j}}<0$, contradicts the linear independence of $\left\{a_{j} \mid j \in A\left(x^{k}\right)\right\}$. Therefore, $\bar{d}^{k} \neq 0$ and as earlier [cf. Eq. (5.10)], we obtain

$$
\nabla f\left(x^{k}\right)^{\prime} \bar{d}^{k} \leq-\frac{1}{2} \bar{d}^{\prime k} \bar{H}^{k} \bar{d}^{k}<0
$$

implying that $\bar{d}^{k}$ is a descent direction. To show that $\bar{d}^{k}$ is also a feasible direction, we must show that

$$
a_{j}^{\prime} \bar{d}^{k} \leq 0, \quad \forall j \in A\left(x^{k}\right)
$$

We know that $a_{j}^{\prime} \bar{d}^{k}=0$ for all $j \in A\left(x^{k}\right)$ except $j=\bar{j}$, and we will show that $a_{j}^{\prime} \bar{d}^{k}<0$. Indeed, by taking inner product of Eq. (5.12) with $\bar{d}^{k}$, we obtain

$$
\nabla f\left(x^{k}\right)^{\prime} \bar{d}^{k}+\sum_{j \in A\left(x^{k}\right)} \mu_{j} a_{j}^{\prime} \bar{d}^{k}=0
$$

and since by construction [cf. Eq. (5.13)], we have $a_{j}^{\prime} \bar{d}^{k}=0$ for all $j \in$ $A\left(x^{k}\right), j \neq \bar{j}$, we obtain

$$
a_{j}^{\prime} \bar{d}^{k}=-\frac{\nabla f\left(x^{k}\right)^{\prime} \bar{d}^{k}}{\mu_{\bar{j}}}
$$

Since $\mu_{\bar{j}}<0$ and $\nabla f\left(x^{k}\right)^{\prime} \bar{d}^{k}<0$, we see that $a_{j}^{\prime} \bar{d}^{k}<0$. Thus, $\bar{d}^{k}$ is a feasible descent direction and a small movement along $\bar{d}^{k}$ makes the $\bar{j}$ th constraint inactive, while it maintains the active status of all other constraints that are active at $x^{k}$.

We have thus completed the description of the manifold suboptimization iteration, together with a method for modifying the manifold of active constraints when no more progress is possible on the current active constraint manifold. It is given by iterations (5.11) or (5.14), depending on whether the scaled projection is done on the manifold of active constraints, or on the manifold obtained after one of the constraints is dropped, respectively.

---

# Positive Definite Quadratic Programming 

Suppose that $f$ is the positive definite quadratic function

$$
f(x)=\frac{1}{2} x^{\prime} Q x+c^{\prime} x
$$

Consider the algorithm of this section with the matrices $H^{k}$ and $\bar{H}^{k}$ in the quadratic subproblems (5.8) and (5.12) chosen to be equal to $Q$. In this case, the quadratic subproblem (5.8) takes the form

$$
\begin{aligned}
& \operatorname{minimize}\left(Q x^{k}+c\right)^{\prime} d+\frac{1}{2} d^{\prime} Q d \\
& \text { subject to } a_{j}^{\prime} d=0, \quad j \in A\left(x^{k}\right)
\end{aligned}
$$

This problem can be viewed as a restricted version of the original quadratic program. Indeed, since $a_{j}^{\prime} x^{k}=b_{j}$ for $j \in A\left(x^{k}\right)$, by letting $x=x^{k}+d$ and by adding the constant term $\frac{1}{2} x^{k^{\prime}} Q x^{k}+c^{\prime} x^{k}$ to its cost function, this problem can be written equivalently as

$$
\begin{aligned}
& \operatorname{minimize} \quad \frac{1}{2} x^{\prime} Q x+c^{\prime} x \\
& \text { subject to } a_{j}^{\prime} x=b_{j}, \quad j \in A\left(x^{k}\right)
\end{aligned}
$$

This is the same problem as the original except that the cost function $f(x)$ is minimized over the currently active manifold $\left\{x \mid a_{j}^{\prime} x=b_{j}, j \in A\left(x^{k}\right)\right\}$ rather than over the full constraint set $\left\{x \mid a_{j}^{\prime} x \leq b_{j}, j=1, \ldots, r\right\}$.

Let now $d^{k}$ be the unique optimal solution of problem (5.16). There are two cases (see Fig. 2.5.4):

1. $d^{k}=0$. Then, one of the active constraints $\bar{j} \in A\left(x^{k}\right)$ involving a negative coefficient $\mu_{\bar{j}}$, must be dropped to form the enlarged manifold of active constraints. Once this is done, we obtain $d^{k}$ as the unique solution of the problem

$$
\begin{aligned}
& \operatorname{minimize}\left(Q x^{k}+c\right)^{\prime} d+\frac{1}{2} d^{\prime} Q d \\
& \text { subject to } a_{j}^{\prime} d=0, \quad j \in A\left(x^{k}\right), j \neq \bar{j}
\end{aligned}
$$

and proceed as in the following case (2).
2. $d^{k} \neq 0$. Then, there are two possibilities:
(a) $x^{k}+d^{k}$ is feasible, in which case we set

$$
x^{k+1}=x^{k}+d^{k}
$$

Then, since problems (5.16) and (5.17) are equivalent, $x^{k+1}$ minimizes $f(x)$ over the current manifold $\left\{x \mid a_{j}^{\prime} x=b_{j}, j \in A\left(x^{k}\right)\right\}$, and in the next iteration, we will drop a constraint as in case (1) above.

---

(b) $x^{k}+d^{k}$ is not feasible, in which case we set

$$
x^{k+1}=x^{k}+\alpha^{k} d^{k}
$$

where

$$
\alpha^{k}=\max \left\{\alpha>0 \mid x^{k}+\alpha d^{k}: \text { feasible }\right\}
$$

Then, at least one (and typically just one) new constraint will be added to the active set for the next iteration.

In a practical algorithm, the operations of adding and dropping constraints from the active set, together with the solution of the corresponding quadratic subproblems (5.17) can be done using efficient linear algebra operations (see [GiM74] and [GMW81] for a more detailed discussion).
![[ch_2_p65_img20.jpeg]]

Figure 2.5.4. Illustration of the quadratic programming algorithm. In (a) the vector $x^{k}+d^{k}$ is feasible and one of the active constraints must be dropped at the next iteration; in (b) the vector $x^{k}+d^{k}$ is not feasible, in which case a new constraint is added to the active set.

We finally show that the quadratic programming algorithm terminates with the unique optimal solution in a finite number of iterations. To this end, let us classify the manifold of active constraints $x^{k}+S\left(x^{k}\right)$ as one of two types: if the unique minimizer of $f$ over $x^{k}+S\left(x^{k}\right)$ is feasible, we say that the manifold is of type $F$ and otherwise we say that the manifold is of type $N$. We note that a manifold of type $F$ can be visited at most once because if $x^{k}$ belongs to such a manifold, then $x^{k+1}$ is the unique minimizer of $f$ over the manifold and a descent algorithm cannot generate the same point twice. On the other hand, each time the algorithm visits a manifold of type $N$, at the next iteration it visits a manifold involving at least one more constraint. Therefore, there can be at most $n$ successive

---

iterations in manifolds of type $N$. Thus, if the number of manifolds of type $F$ is $m$, the number of iterations is bounded by $m n$.

# E XERCISES 

## 5.1

Use the method of this section to solve the three-dimensional quadratic problem

$$
\begin{aligned}
& \operatorname{minimize} f(x)=x_{1}^{2}+2 x_{2}^{2}+3 x_{3}^{2} \\
& \text { subject to } x_{1}+x_{2}+x_{3} \geq 1, \quad 0 \leq x_{1}, 0 \leq x_{2}, 0 \leq x_{3}
\end{aligned}
$$

starting from the point $x^{0}=(0,0,1)$.

## 5.2

Show by example that in the method of this section, if several constraints $a_{j}^{\prime} x \leq b_{j}$ with negative values $\mu_{j}$ are simultaneously dropped, then the corresponding vector $\bar{d}^{k}$ need not be a feasible direction.

### 2.6 AFFINE SCALING FOR LINEAR PROGRAMMING

The gradient projection methods for linear inequality constraints considered so far obtain the descent direction by solving some quadratic programming problem. We saw two extreme types of methods with several others lying in between. In the first type of method, discussed in Sections 2.3 and 2.4 , the projection is done by solving an inequality constrained quadratic problem involving all of the linear constraints; here the set of active constraints may change radically from one iteration to the next. In the second type of method, the manifold suboptimization methods discussed in the preceding section, the projection is done by solving an easier equality constrained quadratic problem; here the set of active constraints typically changes by no more than one constraint per iteration. As a result the method may require a large number of iterations to identify the set of active constraints.

---

In this section we consider a third type of gradient projection method, which has proved particularly interesting for linear programming problems. It shares with the manifold suboptimization method the advantage of solving equality (rather than inequality) constrained quadratic programming problems; however, the optimal solution is approached through the interior of the constraint set, thereby obviating the inefficiency of identifying a large number of active constraints, one-by-one.

We will focus on the linear programming problem

$$
\begin{aligned}
& \operatorname{minimize} c^{\prime} x \\
& \text { subject to } A x=b, \quad x \geq 0
\end{aligned}
$$

where $c \in \Re^{n}$ and $b \in \Re^{m}$ are given vectors, and $A$ is an $m \times n$ matrix of rank $m$. We assume that the optimal cost $\min _{A x=b, x \geq 0} c^{\prime} x$ is finite. Suppose that we have a feasible vector $x^{k}$ such that $x^{k}>0$, that is, $x_{i}^{k}>0$ for $i=1, \ldots, n$ (an initial vector $x_{0}>0$ can be found as discussed in Subsection 2.2.1).

The scaled gradient projection iteration considered in Subsection 2.3.1 [cf. Eqs. (3.12) and (3.13)] is given by

$$
x^{k+1}=x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right)
$$

where $\alpha^{k}$ is a positive stepsize and $\bar{x}^{k}$ solves a quadratic programming problem of the form

$$
\begin{aligned}
& \operatorname{minimize} c^{\prime}\left(x-x^{k}\right)+\frac{1}{2 s^{k}}\left(x-x^{k}\right)^{\prime} H^{k}\left(x-x^{k}\right) \\
& \text { subject to } A x=b, \quad x \geq 0
\end{aligned}
$$

Here, $H^{k}$ is a positive definite matrix and $s^{k}$ is a positive parameter.
We now observe that for a fixed $H^{k}$, the constraints $x \geq 0$ are inactive at the solution $\bar{x}^{k}$ of the quadratic programming problem ( $\bar{x}^{k}>0$ ), provided $s^{k}$ is small enough. Indeed, based on Example 1.5 in Section 2.1, we see that the solution of the equality constrained problem

$$
\begin{aligned}
& \operatorname{minimize} c^{\prime}\left(x-x^{k}\right)+\frac{1}{2 s^{k}}\left(x-x^{k}\right)^{\prime} H^{k}\left(x-x^{k}\right) \\
& \text { subject to } A x=b
\end{aligned}
$$

is given by

$$
\bar{x}^{k}=x^{k}-s^{k}\left(H^{k}\right)^{-1}\left(c-A^{\prime} \lambda^{k}\right)
$$

where

$$
\lambda^{k}=\left(A\left(H^{k}\right)^{-1} A^{\prime}\right)^{-1} A\left(H^{k}\right)^{-1} c
$$

From Eq. (6.4), we see that if $s^{k}$ is sufficiently small, then since $x^{k}>0$, we also have $\bar{x}^{k}>0$, implying that the constraints $x \geq 0$ of the quadratic

---

program (6.2) are superfluous, and that $\bar{x}^{k}$ as given by Eqs. (6.4) and (6.5) solves that program.

By using the expression (6.4) for $\bar{x}^{k}$ in iteration (6.1), and by lumping the scalar $s^{k}$ into the stepsize parameter $\alpha^{k}$, we obtain the iteration

$$
x^{k+1}=x^{k}-\alpha^{k}\left(H^{k}\right)^{-1}\left(c-A^{\prime} \lambda^{k}\right)
$$

where $\lambda^{k}$ is given by Eq. (6.5) and $\alpha^{k}$ is chosen small enough to ensure that $x^{k+1}>0$; see Fig. 2.6.1. In particular we must have $\alpha^{k}<\bar{\alpha}^{k}$ where $\bar{\alpha}^{k}$ is the maximum stepsize for which $x^{k+1}$ is feasible, that is,

$$
\bar{\alpha}^{k}=\max \left\{\alpha \mid x^{k}-\alpha\left(H^{k}\right)^{-1}\left(c-A^{\prime} \lambda^{k}\right) \geq 0\right\}
$$

We know that $\bar{x}^{k}-x^{k}$ is a descent direction, so the cost decrease is proportional to the size of $\alpha^{k}$; this argues for a value of $\alpha^{k}$ close to $\bar{\alpha}^{k}$. On the other hand if $\alpha^{k}$ is very close to $\bar{\alpha}^{k}$, then $x^{k+1}$ will be very close to the constraint boundary, restricting large cost improvements in subsequent iterations. In practice it is common to choose $\alpha^{k}$ between $0.9 \bar{\alpha}^{k}$ and $0.999 \bar{\alpha}^{k}$. Note that $\bar{\alpha}^{k}$ is well-defined as a positive scalar by Eq. (6.7), since if we had $x^{k}-\alpha\left(H^{k}\right)^{-1}\left(c-A^{\prime} \lambda^{k}\right)>0$ for arbitrarily large stepsizes $\alpha$, then we would be able to make the cost $c^{\prime} x$ arbitrarily small, contradicting our assumption that the optimal cost is finite.

An important question has to do with the choice of the scaling matri$\operatorname{ces} H^{k}$. Not every type of scaling matrix can lead to convergence through the interior of the constraint set. For example, it can be seen that choosing $H^{k}$ to be the identity for all $k$, which corresponds to unscaled gradient projection, does not work because for a linear cost, the projected gradient direction will change only if $H^{k}$ changes. What is needed is a matrix with the property that it "bends" the direction $\bar{x}^{k}-x^{k}$ "away from the boundary" of the constraint set.

A particularly interesting choice is

$$
H^{k}=\left(X^{k}\right)^{-2}
$$

where $X^{k}$ is the diagonal matrix having the (positive) coordinates $x_{i}^{k}, i=$ $1, \ldots, n$, along the diagonal. The iteration (6.6) then takes the form

$$
x^{k+1}=x^{k}-\alpha^{k}\left(X^{k}\right)^{2}\left(c-A^{\prime} \lambda^{k}\right)
$$

where

$$
\lambda^{k}=\left(A\left(X^{k}\right)^{2} A^{\prime}\right)^{-1} A\left(X^{k}\right)^{2} c
$$

and is known as the affine scaling method. It was invented in 1967 by Dikin [Dik67], but it received a lot of attention in the West following the proposal by Karmarkar of a different but related method [Kar84]. Karmarkar's method has been shown to have nicer theoretical properties than the affine

---

![[ch_2_p69_img21.jpeg]]

Figure 2.6.1. Illustration of the gradient projection method (6.6) applied to a linear program with three variables and one constraint. The new iterate $x^{k+1}$ is obtained by solving the problem

$$
\begin{aligned}
& \operatorname{minimize} c^{\prime}\left(x-x^{k}\right)+\frac{1}{2 \alpha^{k}}\left(x-x^{k}\right)^{\prime} H^{k}\left(x-x^{k}\right) \\
& \text { subject to } A x=b
\end{aligned}
$$

where $\alpha^{k}$ is such that $x^{k+1}>0$. The scaling matrix $H^{k}$ determines the shape of the ellipsoid corresponding to the iteration. Clearly, $H^{k}$ must change at each iteration, so that the ellipsoid "adapts" to the boundary of the feasible region. The affine scaling method uses the choice $H^{k}=\left(X^{k}\right)^{-2}$, so that if a coordinate $x^{i}$ is small, the corresponding axis of the ellipsoid is also small.
![[ch_2_p69_img22.jpeg]]

Figure 2.6.2. An iteration of the affine scaling method in the original and in the transformed coordinate system under the transformation $y=\left(X^{k}\right)^{-1} x$.

---

scaling method (a polynomial worst-case computational complexity), but in practice the affine scaling method is generally thought to be superior.

The use of the scaling matrix $H^{k}=\left(X^{k}\right)^{-2}$ corresponds to an unscaled gradient projection iteration in the coordinate system $y=\left(X^{k}\right)^{-1} x$. The vector $x^{k}$ is mapped in this system onto the unit vector $y^{k}=(1, \ldots, 1)$. Since the coordinates of $y^{k}$ are uniformly far from the boundary of the constraint $y \geq 0$, this choice of scaling tends to allow large steps and associated cost reductions. Figure 2.6.2 illustrates the iteration in the two coordinate systems.

# Inequality Constrained Linear Programs* 

Consider now the problem

$$
\begin{aligned}
& \text { maximize } b^{\prime} \lambda \\
& \text { subject to } A^{\prime} \lambda \leq c
\end{aligned}
$$

where $\lambda \in \Re^{m}$ is the vector of decision variables, and $b, c$, and $A$ are as in problem (LP1). It will be shown in Section 3.4 that problems (LP1) and (LP2) are dual to each other, in the sense that optimal solutions $\bar{x}$ of (LP1) and $\bar{\lambda}$ of (LP2) are related by $c^{\prime} \bar{x}=b^{\prime} \bar{\lambda}$, but this duality relation will not concern us here. We will develop a version of the affine scaling method, applied to problem (LP2).

Let $w \in \Re^{n}$ be a vector satisfying $A w=b$. [It can be shown that such a vector exists if problem (LP2) has an optimal solution; see Subsection 3.4.2, Prop. 3.4.3, and Example 4.2.] Then by making the transformation $x=c-A^{\prime} \lambda$, we obtain $b^{\prime} \lambda=w^{\prime} A^{\prime} \lambda=w^{\prime}(c-x)$, so problem (LP2) can equivalently be written as

$$
\begin{aligned}
& \operatorname{minimize} w^{\prime} x \\
& \text { subject to } x \geq 0 \text {, and } x=c-A^{\prime} \lambda \text { for some } \lambda \in \Re^{m}
\end{aligned}
$$

Since the constraint $x=c-A^{\prime} \lambda$ for some $\lambda \in \Re^{m}$ defines a linear manifold, this is a linear program of the form (LP1) to which the affine scaling method can be applied. The $k^{\text {th }}$ iteration takes the following form:

Given $x^{k}$ and $\lambda^{k}$ with $x^{k}=c-A^{\prime} \lambda^{k}>0$, we solve the problem

$$
\begin{aligned}
& \operatorname{minimize} w^{\prime}\left(x-x^{k}\right)+\frac{1}{2 s^{k}}\left(x-x^{k}\right)^{\prime} H^{k}\left(x-x^{k}\right) \\
& \text { subject to } x-x^{k}=-A^{\prime}\left(\lambda-\lambda^{k}\right) \text { for some } \lambda \in \Re^{m}
\end{aligned}
$$

where

$$
H^{k}=\left(X^{k}\right)^{-2}
$$

and $X^{k}$ is the diagonal matrix with the coordinates of the vector $x^{k}=$ $c-A^{\prime} \lambda^{k}$ along the diagonal. The affine scaling iteration has the form

$$
x^{k+1}=x^{k}+\alpha^{k}\left(\bar{x}^{k}-x^{k}\right)
$$

---

where $\bar{x}^{k}$ solves problem (6.12) and $\alpha^{k}$ is a positive stepsize that maintains the condition $x^{k+1}>0$. Equivalently, using the relation $A w=b$, problem (6.12) can be written in terms of the vector $\lambda$ as

$$
\begin{aligned}
& \operatorname{minimize}-b^{\prime}\left(\lambda-\lambda^{k}\right)+\frac{1}{2 s^{k}}\left(\lambda-\lambda^{k}\right)^{\prime} A H^{k} A^{\prime}\left(\lambda-\lambda^{k}\right) \\
& \text { subject to } \lambda \in \Re^{m}
\end{aligned}
$$

The solution is

$$
\bar{\lambda}^{k}=\lambda^{k}+s^{k}\left(A H^{k} A^{\prime}\right)^{-1} b
$$

and by using the relation $\bar{x}^{k}-x^{k}=-A^{\prime}\left(\bar{\lambda}^{k}-\lambda^{k}\right)$, the affine scaling iteration (6.13) can be written as

$$
x^{k+1}=x^{k}-\alpha^{k} A^{\prime}\left(A H^{k} A^{\prime}\right)^{-1} b
$$

By using the transformation $x^{k}=c-A^{\prime} \lambda^{k}$, we obtain

$$
A^{\prime} \lambda^{k+1}=A^{\prime} \lambda^{k}+\alpha^{k} A^{\prime}\left(A H^{k} A^{\prime}\right)^{-1} b
$$

By multiplying this iteration by $\left(A A^{\prime}\right)^{-1} A$, we can write it as

$$
\lambda^{k+1}=\lambda^{k}+\alpha^{k}\left(A H^{k} A^{\prime}\right)^{-1} b
$$

where $\alpha^{k}$ is a positive stepsize maintaining the condition $c-A^{\prime} \lambda^{k+1}>0$. This is the affine scaling iteration applied to $\lambda^{k}$.

Because of the duality of problems (LP1) and (LP2) mentioned earlier, iteration (6.16) is called the dual affine scaling method for problem (LP1). In contrast, the iteration (6.8) is called the primal affine scaling method for problem (LP1). However, as we have seen above, the dual method is the same as the primal method applied to the equivalent problem (6.11).

# Convergence Analysis of the Affine Scaling Method* 

The affine scaling method, in the form described above, typically performs very well in practice, but its theoretical convergence to the optimum can be shown only under certain restrictive nondegeneracy assumptions, to be discussed shortly. In particular, let us write the primal affine scaling method (6.8), (6.9) as

$$
x^{k+1}=x^{k}-\beta \bar{\alpha}^{k}\left(X^{k}\right)^{2} z^{k}
$$

where

$$
\begin{gathered}
z^{k}=c-A^{\prime} \lambda^{k} \\
\lambda^{k}=\left(A\left(X^{k}\right)^{2} A^{\prime}\right)^{-1} A\left(X^{k}\right)^{2} c
\end{gathered}
$$

---

$\bar{\alpha}^{k}$ is the maximum stepsize for which the new iterate is feasible [cf. Eq. $(6.7)]$

$$
\begin{aligned}
\bar{\alpha}^{k} & =\max \left\{\alpha \mid x^{k}-\alpha\left(X^{k}\right)^{2} z^{k} \geq 0\right\} \\
& =\max \left\{\left(x_{i}^{k} z_{i}^{k}\right)^{-1} \mid z_{i}^{k}>0, i=1, \ldots, n\right\}
\end{aligned}
$$

and $\beta$ is a scalar in $(0,1)$, which for good performance should be very close to 1 (say $\beta=0.99$ ).

We will derive an expression for the amount of cost improvement at the $k$ th iteration. We have, using the fact $A\left(x^{k}-x^{k+1}\right)=0$,

$$
\begin{aligned}
\Delta \text { cost } & =c^{\prime}\left(x^{k}-x^{k+1}\right) \\
& =\left(c-A^{\prime} \lambda^{k}\right)^{\prime}\left(x^{k}-x^{k+1}\right) \\
& =z^{k^{\prime}}\left(\beta \bar{\alpha}^{k}\left(X^{k}\right)^{2} z^{k}\right) \\
& =\beta \bar{\alpha}^{k}\left\|X^{k} z^{k}\right\|^{2}
\end{aligned}
$$

and by using the expression (6.20) for $\bar{\alpha}^{k}$, we have

$$
\Delta \text { cost }=\beta\left\|X^{k} z^{k}\right\|^{2} \max \left\{\left(x_{i}^{k} z_{i}^{k}\right)^{-1} \mid z_{i}^{k}>0, i=1, \ldots, n\right\}
$$

Combining this relation with the fact

$$
\begin{aligned}
\left\|X^{k} z^{k}\right\| & =\left(\sum_{i=1}^{n}\left(x_{i}^{k} z_{i}^{k}\right)^{2}\right)^{1 / 2} \\
& \geq \min \left\{x_{i}^{k} z_{i}^{k} \mid z_{i}^{k}>0, i=1, \ldots, n\right\} \\
& =\frac{1}{\max \left\{\left(x_{i}^{k} z_{i}^{k}\right)^{-1} \mid z_{i}^{k}>0, i=1, \ldots, n\right\}}
\end{aligned}
$$

we finally obtain

$$
\begin{aligned}
\Delta \text { cost } & \geq \beta\left\|X^{k} z^{k}\right\| \\
& =\beta\left(\sum_{i=1}^{n}\left(x_{i}^{k} z_{i}^{k}\right)^{2}\right)^{1 / 2} \\
& \geq \frac{\beta}{2} \sum_{i=1}^{n}\left|x_{i}^{k} z_{i}^{k}\right|
\end{aligned}
$$

Thus the cost improvement at the $k$ th iteration is at least proportional to each of the products $\left|x_{i}^{k} z_{i}^{k}\right|$, implying that

$$
x_{i}^{k} z_{i}^{k} \rightarrow 0, \quad i=1, \ldots, n
$$

Suppose that the sequence $\left\{\left(x^{k}, \lambda^{k}\right)\right\}$ converges to some $(\bar{x}, \bar{\lambda})$, and let $\bar{z}=c-A^{\prime} \bar{\lambda}$. Then from Eq. (6.23), we have

$$
\bar{x}_{i} \bar{z}_{i}=0, \quad i=1, \ldots, n
$$

---

from the feasibility of $x^{k}$, we have

$$
A \bar{x}=b, \quad \bar{x} \geq 0
$$

and from Eq. (6.17), we have

$$
\bar{z} \geq 0
$$

[if $\bar{z}_{1}<0$ for some $i$, then from Eq. (6.17) it would follow that $x_{i}^{k+1}>x_{i}^{k}>$ 0 for all $k$, implying $\bar{x}_{i}>0$ and contradicting Eq. (6.24)].

It can be shown that Eqs. (6.24)-(6.26) imply that $\bar{x}$ is an optimal solution of the primal problem (LP1), and $\bar{\lambda}$ is an optimal solution of the dual problem (LP2). This is a consequence of the duality theory to be developed in Subsection 3.4.2 (particularly Example 4.2); we refer the reader to that section since the optimality of $\bar{x}$ and $\bar{\lambda}$ will not be used further in this section.

We have thus shown that if the sequence $\left\{\left(x^{k}, \lambda^{k}\right)\right\}$ converges to some $(\bar{x}, \bar{\lambda})$, then $\bar{x}$ and $\bar{\lambda}$ are optimal primal and dual solutions, respectively. Unfortunately, however, it is not possible to guarantee the convergence of $\left\{\left(x^{k}, \lambda^{k}\right)\right\}$ without additional quite restrictive assumptions generally known as nondegeneracy. One such assumption is that for each primal feasible solution $x$ we must have $x_{i}>0$ for at least $m$ indices $i$, and that for each dual feasible solution $\lambda$ and corresponding vector $z=c-A^{\prime} \lambda$, we must have $z_{i}>0$ for at least $n-m$ indices $i$. Another type of nondegeneracy assumption under which convergence can be proved is developed in Exercise 6.2. Such nondegeneracy assumptions are seldom satisfied in practice, but have often been used to simplify various linear programming analyses. If we do not assume nondegeneracy, it is necessary to modify the stepsize procedure of the affine scaling method in order to prove convergence to the optimum (see the references). Unfortunately, while the theoretical convergence properties of the method are then improved, the practical performance becomes substantially slower.

# E XERCISES 

## 6.1

Consider the linear program

$$
\begin{aligned}
& \text { minimize } x_{1}+2 x_{2}+3 x_{3} \\
& \text { subject to } x_{1}+x_{2}+x_{3}=1, \quad x \geq 0
\end{aligned}
$$

---

Write a computer program to implement the affine scaling method for this problem. Solve the problem for the starting points $x^{0}=(.8, .15, .05)$ and $x^{0}=(.1, .2, .7)$.

# 6.2 (Convergence of the Affine Scaling Method) 

This exercise deals with the convergence of the primal affine scaling method (6.8), (6.9) under a nondegeneracy assumption, and requires some knowledge of the theory of polyhedral convex sets given in Sections B. 3 and B. 4 of Appendix B. Show that $\left\{x^{k}\right\}$ converges to the optimal primal solution $x^{*}$ and $\left\{\lambda^{k}\right\}$ converges to an optimal dual solution $\lambda^{*}$, assuming the following:
(1) Every vector $z$ of the form $z=c-A^{\prime} \lambda$, where $\lambda \in \Re^{m}$, has at most m zero coordinates.
(2) The extreme points of the polyhedron $X=\{x \mid A x=b, x \geq 0\}$ have distinct costs and the primal optimal solution is unique.
(3) The matrix $A X^{2} A^{\prime}$ is nonsingular for all $x \in X$.

Hint: Complete the details in the following argument. Assumption (3) implies that the inverse $\left(A X^{2} A^{\prime}\right)^{-1}$ is continuous over $X$, and that if $x \in X$, then $x$ must have at most $n-m$ nonzero coordinates. Assumption (2) implies that for every scalar $\gamma$, the level set $\{x \mid c^{\prime} x \leq \gamma, A x=b, x \geq 0\}$ is bounded, so that $\left\{x^{k}\right\}$ is bounded, and therefore also that $\left\{\lambda^{k}\right\}$ and $\left\{z^{k}\right\}$ (where $z^{k}=$ $c-A^{\prime} \lambda^{k}$ ) are bounded. If $(x, z)$ is any limit point of $\left\{\left(x^{k}, z^{k}\right)\right\}$, we have $X z=0$. Then $x$ has exactly $n-m$ nonzero coordinates, and $z$ has exactly $m$ nonzero coordinates, implying that $x$ is an extreme point of the polyhedron $X$. Using assumption (2), conclude that $\left\{x^{k}\right\}$ converges, implying that $\left\{\lambda^{k}\right\}$ also converges.

### 2.7 BLOCK COORDINATE DESCENT METHODS*

We briefly discussed coordinate descent methods for unconstrained minimization in Section 1.7. These methods can be generalized to solve the problem

$$
\begin{aligned}
& \operatorname{minimize} f(x) \\
& \text { subject to } x \in X
\end{aligned}
$$

where $X$ is a Cartesian product of closed, convex sets $X_{1}, \ldots, X_{m}$ :

$$
X=X_{1} \times X_{2} \times \cdots \times X_{m}
$$

We assume that $X_{i}$ is a closed convex subset of $\Re^{n_{i}}$ and $n=n_{1}+\cdots+n_{m}$. The vector $x$ is partitioned as

$$
x=\left(x_{1}, x_{2}, \ldots, x_{m}\right)
$$

---

where each subvector $x_{i}$ belongs to $\Re^{n_{i}}$, so the constraint $x \in X$ is equivalent to

$$
x_{i} \in X_{i}, \quad i=1, \ldots, m
$$

Let us assume that for every $x \in X$ and every $i=1, \ldots, m$, the optimization problem

$$
\begin{aligned}
& \operatorname{minimize} f\left(x_{1}, \ldots, x_{i-1}, \xi, x_{i+1}, \ldots, x_{m}\right) \\
& \text { subject to } \xi \in X_{i}
\end{aligned}
$$

has at least one solution. The following algorithm, known as block coordinate descent or nonlinear Gauss-Seidel method, generates the next iterate $x^{k+1}=\left(x_{1}^{k+1}, \ldots, x_{m}^{k+1}\right)$, given the current iterate $x^{k}=\left(x_{1}^{k}, \ldots, x_{m}^{k}\right)$, according to the iteration

$$
x_{i}^{k+1}=\arg \min _{\xi \in X_{i}} f\left(x_{1}^{k+1}, \ldots, x_{i-1}^{k+1}, \xi, x_{i+1}^{k}, \ldots, x_{m}^{k}\right), \quad i=1, \ldots, m
$$

Thus, at each iteration, the cost is minimized with respect to each of the "block coordinate" vectors $x_{i}^{k}$, taken in cyclic order. Naturally, the method makes practical sense if the minimization in Eq. (7.3) is fairly easy. This is frequently so when each $x_{i}$ is a scalar, but there are also other cases of interest, where $x_{i}$ is a multidimensional vector.

The following proposition gives the basic convergence result for the method. Note that there is a strict convexity assumption in this proposition, which is not required to prove corresponding convergence results for gradient methods (cf. Prop. 1.2.1 in Section 1.2) and for feasible direction methods (cf. Prop. 2.2.1 in Section 2.2). The convexity assumption on $f$ is not absolutely essential in order to assert convergence of the block coordinate descent method. It turns out, however, that it is necessary to make an assumption implying that the minimum in Eq. (7.3) is uniquely attained (see [Zan69'). The need for this is not obvious but has been demonstrated by an example given in [Pow73].

Proposition 2.7.1: (Convergence of Block Coordinate Descent) Suppose that $f$ is continuously differentiable and convex over the set $X$ of Eq. (7.1). Furthermore, suppose that for each $i, f$ is a strictly convex function of $x_{i}$, when the values of the other components of $x$ are held constant. Let $\left\{x^{k}\right\}$ be the sequence generated by the block coordinate descent method (7.3). Then, every limit point of $\left\{x^{k}\right\}$ minimizes $f$ over $X$.

Proof: Let

$$
z_{i}^{k}=\left(x_{1}^{k+1}, \ldots, x_{i}^{k+1}, x_{i+1}^{k}, \ldots, x_{m}^{k}\right)
$$

---

Using the definition (7.3) of the method, we obtain

$$
f\left(x^{k}\right) \geq f\left(z_{1}^{k}\right) \geq f\left(z_{2}^{k}\right) \geq \cdots \geq f\left(z_{m-1}^{k}\right) \geq f\left(x^{k+1}\right), \quad \forall k
$$

Let $\bar{x}=\left(\bar{x}_{1}, \ldots, \bar{x}_{m}\right)$ be a limit point of the sequence $\left\{x^{k}\right\}$. Notice that $\bar{x} \in X$ because $X$ is closed. Equation (7.4) implies that the sequence $\left\{f\left(x^{k}\right)\right\}$ converges to $f(\bar{x})$. It now remains to show that $\bar{x}$ minimizes $f$ over $X$.

Let $\left\{x^{k_{j}} \mid j=0,1, \ldots\right\}$ be a subsequence of $\left\{x^{k}\right\}$ that converges to $\bar{x}$. We first show that $\left\{x_{1}^{k_{j}+1}-x_{1}^{k_{j}}\right\}$ converges to zero as $j \rightarrow \infty$. Assume the contrary or, equivalently, that $\left\{z_{1}^{k_{j}}-x^{k_{j}}\right\}$ does not converge to zero. Let $\gamma^{k_{j}}=\left\|z_{1}^{k_{j}}-x^{k_{j}}\right\|$. By possibly restricting to a subsequence of $\left\{k_{j}\right\}$, we may assume that there exists some $\bar{\gamma}>0$ such that $\gamma^{k_{j}} \geq \bar{\gamma}$ for all $j$. Let $s_{1}^{k_{j}}=\left(z_{1}^{k_{j}}-x^{k_{j}}\right) / \gamma^{k_{j}}$. Thus, $z_{1}^{k_{j}}=x^{k_{j}}+\gamma^{k_{j}} s_{1}^{k_{j}},\left\|s_{1}^{k_{j}}\right\|=1$, and $s_{1}^{k_{j}}$ differs from zero only along the first block-component. Notice that $s_{1}^{k_{j}}$ belongs to a compact set and therefore has a limit point $\bar{s}_{1}$. By restricting to a further subsequence of $\left\{k_{j}\right\}$, we assume that $s_{1}^{k_{j}}$ converges to $\bar{s}_{1}$.

Let us fix some $\epsilon \in[0,1]$. Notice that $0 \leq \epsilon \bar{\gamma} \leq \gamma^{k_{j}}$. Therefore, $x^{k_{j}}+\epsilon \bar{\gamma} s_{1}^{k_{j}}$ lies on the segment joining $x^{k_{j}}$ and $x^{k_{j}}+\gamma^{k_{j}} s_{1}^{k_{j}}=z_{1}^{k_{j}}$, and belongs to $X$ because $X$ is convex. Using the convexity of $f$, and the fact that $z_{1}^{k_{j}}$ minimizes $f$ over all $x$ that differ from $x^{k_{j}}$ along the first block-component, we obtain

$$
f\left(z_{1}^{k_{j}}\right)=f\left(x^{k_{j}}+\gamma^{k_{j}} s_{1}^{k_{j}}\right) \leq f\left(x^{k_{j}}+\epsilon \bar{\gamma} s_{1}^{k_{j}}\right) \leq f\left(x^{k_{j}}\right)
$$

Since $f\left(x^{k}\right)$ converges to $f(\bar{x})$, Eq. (7.4) shows that $f\left(z_{1}^{k}\right)$ also converges to $f(\bar{x})$. We now take the limit as $j$ tends to infinity, to obtain $f(\bar{x}) \leq$ $f\left(\bar{x}+\epsilon \bar{\gamma} \bar{s}_{1}\right) \leq f(\bar{x})$. We conclude that $f(\bar{x})=f\left(\bar{x}+\epsilon \bar{\gamma} \bar{s}_{1}\right)$, for every $\epsilon \in[0,1]$. Since $\bar{\gamma} \bar{s}_{1} \neq 0$, this contradicts the strict convexity of $f$ as a function of the first block-component. This contradiction establishes that $x_{1}^{k_{j}+1}-x_{1}^{k_{j}}$ converges to zero. In particular, $z_{1}^{k_{j}}$ converges to $\bar{x}$.

From the definition (7.3) of the algorithm, we have

$$
f\left(z_{1}^{k_{j}}\right) \leq f\left(x_{1}, x_{2}^{k_{j}}, \ldots, x_{m}^{k_{j}}\right), \quad \forall x_{1} \in X_{1}
$$

Taking the limit as $j$ tends to infinity, we obtain

$$
f(\bar{x}) \leq f\left(x_{1}, \bar{x}_{2}, \ldots, \bar{x}_{m}\right), \quad \forall x_{1} \in X_{1}
$$

Using the conditions for optimality over a convex set (Prop. 2.1.2 in Section 2.1), we conclude that

$$
\nabla_{1} f(\bar{x})^{\prime}\left(x_{1}-\bar{x}_{1}\right) \geq 0, \quad \forall x_{1} \in X_{1}
$$

---

where $\nabla_{i} f$ denotes the gradient of $f$ with respect to the component $x_{i}$.
Let us now consider the sequence $\left\{z_{1}^{k_{j}}\right\}$. We have already shown that $z_{1}^{k_{j}}$ converges to $\bar{x}$. A verbatim repetition of the preceding argument shows that $x_{2}^{k_{j}+1}-x_{2}^{k_{j}}$ converges to zero and $\nabla_{2} f(\bar{x})^{\prime}\left(x_{2}-\bar{x}_{2}\right) \geq 0$ for every $x_{2} \in X_{2}$. Continuing inductively, we obtain $\nabla_{i} f(x)^{\prime}\left(x_{i}-\bar{x}_{i}\right) \geq 0$ for every $x_{i} \in X_{i}$ and for every $i$. Adding these inequalities, and using the Cartesian product structure of the set $X$, we conclude that $\nabla f(\bar{x})^{\prime}(x-\bar{x}) \geq 0$ for every $x \in X$. In view of the convexity of $f$, this shows that $\bar{x}$ minimizes $f$ over the set $X$ (Prop. 2.1.2 in Section 2.1). Q.E.D.

# E XERCISES 

## 7.1 (The Proximal Minimization Algorithm)

Let $f: \Re^{n} \mapsto \Re$ be a continuously differentiable convex function, let $X$ be a closed convex set, and let $c$ be a positive scalar.
(a) Show that the algorithm (cf. Exercise 3.8 in Section 2.3 and Subsection 5.4.6)

$$
x^{k+1}=\arg \min _{x \in X}\left\{f(x)+\frac{1}{2 c}\left\|x-x^{k}\right\|^{2}\right\}
$$

is a special case of the block coordinate descent method applied to the problem

$$
\begin{aligned}
& \text { minimize } f(x)+\frac{1}{2 c}\|x-y\|^{2} \\
& \text { subject to } x \in X, \quad y \in \Re^{n}
\end{aligned}
$$

which is equivalent to the problem of minimizing $f$ over $X$.
(b) Derive a convergence result based on Prop. 2.7.1 for the algorithm of part (a).
(c) Show that if $f$ has at least one minimizing point over $X$, the entire sequence $\left\{x^{k}\right\}$ converges to some such point. Hint: We have by definition

$$
f\left(x^{k+1}\right)+\frac{1}{2 c}\left\|x^{k+1}-x^{k}\right\|^{2} \leq f(x)+\frac{1}{2 c}\left\|x-x^{k}\right\|^{2}, \quad \forall x \in X
$$

so that $\left\|x^{k+1}-x^{k}\right\| \leq\left\|x-x^{k}\right\|$ for all $x$ in the set

$$
X^{k}=\left\{x \in X \mid f(x) \leq f\left(x^{k+1}\right)\right\}
$$

Hence $x^{k+1}$ is the unique projection of $x^{k}$ on $X^{k}$, and we have $\left(x^{k+1}-\right.$ $\left.x^{k}\right)^{\prime}\left(x-x^{k}\right) \geq 0$ for all $x \in X^{k}$. Conclude that for every $x^{*}$ that minimizes $f$ over $X$, we have $\left\|x^{*}-x^{k+1}\right\| \leq\left\|x^{*}-x^{k}\right\|$.

---

# 7.2 (Parallel Projections Algorithm) 

We are given $m$ closed convex sets $C_{1}, C_{2}, \ldots, C_{m}$ in $\Re_{n}$, and we want to find a point in their intersection. Consider the equivalent problem

$$
\begin{array}{ll}
\operatorname{minimize} & \frac{1}{2} \sum_{i=1}^{m}\left\|x_{i}-x\right\|^{2} \\
\text { subject to } & x \in \Re^{n}, x_{i} \in C_{i}, i=1, \ldots, m
\end{array}
$$

where the variables of the optimization are $x, x_{1}, \ldots, x_{m}$. Derive a block coordinate descent algorithm involving projections on each of the sets $C_{i}$ that can be carried out independently of each other. State a convergence result for this algorithm.

### 2.8 NOTES AND SOURCES

Section 2.2: Feasible direction methods were first systematically investigated in [Zou60]. Convergence analyses and additional methods may be found in [ToV67], [Pol71], [PiP73], and [Zou76]. The conditional gradient method was first proposed for quadratic programs in [FrW56]. Extensive discussions of its application to more general problems is found in [LeP65] and [DeR70]. The convergence rate of the method was analyzed in [CaC68], [Dun79], [Dun80a], and [DuS83].

Section 2.3: Gradient projection methods with a constant stepsize were independently proposed in [Gol64] and [LeP65]. The Armijo rule along the projection arc was proposed in [Ber76c]. For convergence analysis, see [Dun81], [GaB82], [GaB84], [CaM87], [Dun87], and [Dun88a]. For analysis of various aspects of the convergence rate, see [Ber76c], [Dun81], [BeG82], [DuS83], [Dun87], [Tse91a], [LuT92b], [LuT93a], [LuT93b], [LuT94b], and [GiK95]. Surveys are given in [Dun88a], [Dun94]. Variations of gradient projection methods that are based on constraint identification are given in [Ber76c], [DeT83], [MoT89], [BMT90], and [MMZ95]. The gradient projection method is particularly well-suited for large-scale problems with relatively simple constraint structure. Examples of such problems arise in optimal control (see e.g., [Pol73], [Ber76c], [Ber82b], [Dun88a], [Dun91a], [Dun94]), and in multicommodity flow problems from communications and transportation (see e.g., [BeG83], [BeT89], [BeG92], [LuT94b]). Newton's method for constrained problems is discussed in [Dun80b] and [HuD84].
Section 2.4: Two metric-projection methods were proposed in [Ber82b] and were generalized in [GaB84]; see also [BeG83]. Subsequent works include [Dun88b], [GaD88], [Dun91a], [Dun91b], [Dun93a], [Dun93b], and [Dun94].

---

Section 2.5: Manifold suboptimization methods draw their origin from the gradient projection method of Rosen [Ros60b]. Extensive discussions of these methods can be found in [GiM74] and [GMW81]. A particular method of this type, the reduced gradient method, has been used as the basis for general purpose nonlinear programming software [LaW78]. A recent convergence analysis of Rosen's gradient projection method is given in [DuZ89]. The application of manifold suboptimization methods to quadratic programming is discussed in [Zan69]. Efficient implementations are given in [GiM74], [GMW81], and [GoI82].

Manifold suboptimization methods relate to one of the oldest and most popular optimization algorithms, the simplex method for linear programming. This method applies to problems of the form

$$
\begin{aligned}
& \text { minimize } c^{\prime} x \\
& \text { subject to } x \in X
\end{aligned}
$$

where $X$ is a polyhedron in $\Re^{n}$ that has at least one vertex (extreme point) and $c$ is a given vector. The idea of the simplex method is to start at some vertex of $X$ and then to generate a sequence of vertices of $X$. These vertices satisfy two properties:
(a) Two successive vertices are connected by an edge; that is, they both lie on some one-dimensional active constraint manifold.
(a) Each vertex has a lower cost than the preceding vertices.

Assuming the linear cost function attains a minimum over $X$, the method is guaranteed to terminate in a finite number of iterations at an optimal vertex. The reason is that $X$ has a finite number of vertices (as shown in Appendix B), a vertex cannot be repeated because of property (b) above, and according to the fundamental theorem of linear programming [Prop. B. 21 (c) in Appendix B], if the linear cost function attains a minimum over $X$, then this minimum is attained at some vertex of $X$.

The selection of the next vertex starting from a given vertex in the simplex method is guided by iterative descent; that is, out of all the edges that lead to neighboring vertices, one that corresponds to a descent direction is chosen. If we view a vertex as an active constraint manifold of dimension 0 and an edge as an active constraint manifold of dimension 1, it is seen that the movement to a neighboring vertex of lower cost encapsulates the basic operations of the manifold suboptimization method: first dropping a constraint at a vertex, then moving along a descent direction on the corresponding one-dimensional manifold (edge), until a new constraint is encountered and the active manifold becomes a vertex again. The choice of a constraint to drop is identical with the one of the manifold suboptimization method (except for streamlining the linear algebra). Thus the behavior of the manifold suboptimization method, when it is applied to linear programming problems and it is started at a vertex of the constraint

---

polyhedron, is identical to the one of the simplex method. The versions of the simplex method used in practice involve a lot of important implementation technology, which is beyond the scope of this book (see, e.g., [GMW91]). This technology is needed to enhance efficiency and also to deal with degeneracy (more than $n$ linear constraints active at a given vertex). Manifold suboptimization methods can benefit from the use of some of this technology. We refer to the literature (e.g., [GiM74], [GMW81]) for more details.

Section 2.6: The affine scaling method was first proposed in [Dik67], and was rediscovered many years later in [Bar86] and [VMF86], following the emergence of interior point methods as serious competitors to the simplex method for linear programming applications. The convergence of the method without restrictive nondegeneracy assumptions is analyzed in [Tsu91], [TsL92], [MTW93], and [TsM95]. Note that while we have focused on the linear programming case, the ideas of the affine scaling method apply to problems with a more general cost function.
Section 2.7: The convergence to a unique limit and the rate of convergence of the coordinate descent method are analyzed in [LuT91] and [LuT92a]. Coordinate descent methods are often well suited for solving dual problems (see Section 6.2), and within this specialized context there has been much convergence analysis (see the references given for Section 6.2).