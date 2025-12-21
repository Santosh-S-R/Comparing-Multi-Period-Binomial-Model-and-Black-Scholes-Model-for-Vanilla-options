# Motivation

After suitable condensation, the price of a vanilla option under the multi-period binomial model can be expressed using the **Cox–Ross–Rubinstein (CRR) formula**, shown below:

<img src="images/crr_exp.png" alt="Cox–Ross–Rubinstein equation" width="450">

One may immediately notice the strong similarity between this expression and the **Black–Scholes formula**:

<img src="images/BS_eq.png" alt="Black–Scholes equation" width="450">

This work shows the convergence of the models under certain conditions, why this happens, and explores some interesting background, applications, and optimisations.

---

## Project Objectives

**Legend**  
_italic_: Joint contribution  
**Bold**: Individual contribution

My team and I structured our work around the following goals:

- Introduce the core financial instruments, with a detailed exposition of the structure of European vanilla options.
- _Present the key mathematical foundations underlying both the Black–Scholes model and the binomial model_.
- _Show the binomial model option price at an granular level, explicitly demonstrating its convergence to the analytical Black–Scholes solution_.
- **Implement the Black–Scholes model in C++ and numerically demonstrate the convergence of the Cox–Ross–Rubinstein model as the number of time steps increases
  This includes a numerical optimisation to prevent factorial values from reaching type limits.**
- **Identify, analyse, and present key observations arising from the numerical experiments.**
- **Conduct a literature review, including a discussion of the trinomial model.**
- Explore extensions to more complex instruments, specifically American options.
