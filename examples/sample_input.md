# Chapter 3: Linear Algebra

## 3.1 Vector Spaces

A **vector space** $V$ over a field $F$ is a set equipped with two operations: vector addition and scalar multiplication, satisfying the following axioms:

### Definition 3.1.1
Let $V$ be a non-empty set and $F$ be a field. Then $V$ is called a vector space over $F$ if:

1. **Closure under addition**: For all $\mathbf{u}, \mathbf{v} \in V$, we have $\mathbf{u} + \mathbf{v} \in V$
2. **Associativity**: $(\mathbf{u} + \mathbf{v}) + \mathbf{w} = \mathbf{u} + (\mathbf{v} + \mathbf{w})$ for all $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$
3. **Additive identity**: There exists $\mathbf{0} \in V$ such that $\mathbf{v} + \mathbf{0} = \mathbf{v}$ for all $\mathbf{v} \in V$
4. **Additive inverse**: For each $\mathbf{v} \in V$, there exists $-\mathbf{v} \in V$ such that $\mathbf{v} + (-\mathbf{v}) = \mathbf{0}$

### Theorem 3.1.2
Let $V$ be a vector space over field $F$. Then:
$$\mathbf{0} \cdot \mathbf{v} = \mathbf{0}$$
for all $\mathbf{v} \in V$, where $\mathbf{0}$ on the left is the zero scalar and $\mathbf{0}$ on the right is the zero vector.

**Proof**: Let $\mathbf{v} \in V$. Then:
$$\mathbf{0} \cdot \mathbf{v} = (\mathbf{0} + \mathbf{0}) \cdot \mathbf{v} = \mathbf{0} \cdot \mathbf{v} + \mathbf{0} \cdot \mathbf{v}$$

Adding the additive inverse of $\mathbf{0} \cdot \mathbf{v}$ to both sides:
$$\mathbf{0} = \mathbf{0} \cdot \mathbf{v}$$

## 3.2 Linear Independence

### Definition 3.2.1
A set of vectors $\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n\}$ in vector space $V$ is called **linearly independent** if the only solution to:
$$c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_n\mathbf{v}_n = \mathbf{0}$$
is $c_1 = c_2 = \cdots = c_n = 0$.

### Example 3.2.2
Consider the vectors in $\mathbb{R}^3$:
- $\mathbf{v}_1 = (1, 0, 0)$
- $\mathbf{v}_2 = (0, 1, 0)$ 
- $\mathbf{v}_3 = (0, 0, 1)$

These vectors are linearly independent because if:
$$c_1(1, 0, 0) + c_2(0, 1, 0) + c_3(0, 0, 1) = (0, 0, 0)$$

Then $(c_1, c_2, c_3) = (0, 0, 0)$, which implies $c_1 = c_2 = c_3 = 0$.

## 3.3 Basis and Dimension

### Definition 3.3.1
A **basis** for a vector space $V$ is a set of vectors that is both linearly independent and spans $V$.

### Theorem 3.3.2 (Basis Extension Theorem)
Every linearly independent set in a finite-dimensional vector space can be extended to a basis.

**Proof Sketch**: Let $S = \{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ be linearly independent in $V$. If $S$ spans $V$, we're done. Otherwise, there exists $\mathbf{w} \in V$ not in $\text{span}(S)$. Then $S \cup \{\mathbf{w}\}$ is linearly independent. Continue this process until we obtain a spanning set.
