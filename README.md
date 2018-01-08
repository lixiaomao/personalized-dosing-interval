# Determining Personalized Dosing Interval


This package for determining Personalized Dosing Interval in the cases where there is one continuous treatment but the dosages needs to personalized. 
we propose two approaches, Probability Dosing Interval(PDI) and Expectation Dosing Interval(EDI). The former yields a potential outcome which is greater than a specific value with a certain probability, the latter allows the expectation of the potential outcome to be greater than a specific value. 
According to the assumption of the potential outcome, one can choose to estimate the lower boundary, upper boundary, or both boundary of the interval. 

There are 3 methods to choose: 

- DC-based algorithms use Difference-in-Difference algorithm and solves the non-convex problem by iteratively solving a convex relaxation. There are three fitters, penalized quantile regression(rq), support vector machine(svmLinear and svmRadial),  and reinforcement learning tree(RLT). Global solution is not guaranteed for these approaches.
- Tree-based algorithm uses no relaxation and solves the problem directly. Global solution is guaranteed to a certain precision.  
- Ordinal approach takes the continuous treatment as ordinal treatments with several levels. 

This work was presented in JSM 2017. 