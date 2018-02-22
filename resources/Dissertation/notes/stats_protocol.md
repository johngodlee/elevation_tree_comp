A protocol for data exploration to avoid common statistical problems

Zuur et al. 2010

-   Provides a protocol for data exploration, discusses current tools to
    detect outliers, heterogeneity of variances, collinearity,
    dependence of observations, problems with interactions, double
    zeroes in multivariate analysis, zero inflation in generalize linear
    modelling, correct types of relationship between dependent and
    independent variables

-   ­­It is important to remember that the 'linear' regression does not
    refer to the linearity of the relationship between x and y, but to
    the parameters used in the model. ‘linear’ regressions can easily
    model non-linear relationships using quadratics and interactions

-   Decisions about what models to test should be made before everything
    else, based on biological understanding of the system (Burnham &
    Anderson, 2002)

-   They recommend this workflow

    1.  Identify Outliers – Box plot and Cleveland dotplot

    2.  Homogeneity of y – Conditional boxplot

    3.  Normality of y – Histogram or Q-Q plot

    4.  Zero troubles in y – Frequency plot or corrgram

    5.  Collinearity of x – VIF & scatterplots & PCA

    6.  Relationships of y & x – multi-panel scatterplots and
        conditional boxplots

    7.  Interactions – coplots

    8.  Intependence of y – ACF & variogram, plot y against time/space

-   Laara 2009 gives seven reasons for not applying preliminary tests
    for normality

1.  Outliers in x and y

    a.  Outliers are treated differently by different statistical
        frameworks

    b.  Boxplots, but Cleveland dotplots are better at identifying
        ACTUAL outliers

    c.  In regression, outliers in the response are harder to deal with,
        better to choose a probability distribution that allows greater
        variation for large mean values, e.g. gamma for continuous data,
        poisson or negative binomial for count data, because it
        preserves the original data

2.  Homogeneity of variances in y (/homoscedasticity in
    linear regression)

    a.  Important in regression and ANOVA

    b.  Should assess in regression type models by plotting the
        residuals vs the fitted values, the residual variation should be
        similar across the fitted range

    c.  If homogeneity of residuals cannot be attained you could use a
        generalized least squares model (Pinheiro & Bates 2000, Zuur
        2009a

3.  Normal distribution?

    a.  If you want to perform a model to assess whether there is
        significant sepearation between groups, you must assess
        normality of the variable *within each group*!

    b.  In linear regression we assume normality of all observations at
        a particular x value – but unless we have lots (&gt;25)
        observations per x value, we can’t do this, instead we use a
        histogram of pooled residuals across the x range

        i.  It is important to check that normality is hed
            within groups. Apparent skewness can actually just be caused
            by mean differences between groups, in this case
            transformation is not advisable as it only reduces the
            ability to delineate differences between groups

4.  Zeroes!

    a.  GLM analysis is appropro

    b.  Zero-inflated GLMs if there are lots of zeroes (zuur et al.
        2009a, Cameron and Trivedi 1998)

5.  Collinearity among x’s

    a.  Can test by making a table of all predictors and their p-values
        & VIFs

    b.  $VIF = 1/(1 - R_{j}^{2})$, where R^2^ = the R^2^ from a linear
        regression model in which covariate X~j~ is used as a response
        variable, and all other covariates as explanatory
        variables, i.e. it shows how much variation in X~j~ is explained
        by all other covariates. High VIF = high collinearity

    c.  Sequentially drop the covariates with the highest VIF, until the
        VIFs of all others drop below 10, 3, etc.

    d.  Collinearity is especially important when ecological signals are
        weak

    e.  Could also do:

        i.  PCA biplot (Jolliffe 2002) of all covariates

        ii. Correlation coefficients

6.  Relationships between x and y

    a.  Important to remember that just because the scatterplots have no
        clear patterns of change in y with x, it doesn’t mean that there
        isn’t a significant relationship between y and x+\*z
        (multiple predictors)

7.  –

8.  Are observations of response (y) independent

    a.  Could account for this in the modelling structure

        i.  Using spatial or temporal variables as a predictor
            (explicitly model them)

        ii. Nesting data in a hierarchical structure

        iii. Random effects in mixed effects model framework (Pinheiro &
            Bates 2000)

    b.  Can plot the response against spatial and temporal covariates

        i.  If there is no clear sequence of observations you can apply
            a dependence structure using random effects

    c.  Can plot an autocorrelation function to show what the lag of
        dependence is

    d.  e.  


