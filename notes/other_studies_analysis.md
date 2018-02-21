How have other studies across elevational gradients done their analysis?

Asner et al. 2014

-   Used means and averages mostly

    -   Not trying to answer a question so know need for hypothesis
        testing

-   Used Pearson’s Correlation

-   R\^2 of how much elevation accounted for differences in canopy
    variables

-   Linear regression of elevation causing changes in canopy variables

Ettinger *et al.* 2013

-   This study seems fairly similar to mine

-   Elevational gradient but across Mt Rainier in Washington

-   Linear Mixed effects model – to evaluate relationship between growth
    and climate for each species and life stage at each sampling
    location

    -   Allows differences among individual tree responses to climate to
        be accommodated

    -   Individual tree and year as random effects

        -   To account for non-independence of data from the same
            individual or within years (Crawley 2007)

        -   Random slopes structure for individual tree

        -   Intercept only structure for year (Zuur et al. 2009)

    -   All climate variables were fixed effects

    -   One model per species

    -   Used AIC to compare the best fitting model

        -   use a null model against the best model to see how much
            variation the climate explains in growth
            (AIC~null~-AIC~Best~)

    -   Calculated the significance of coefficients using “LanguageR”
        package, estimates P values usinhg MCMC sampling (Baayen et al.
        2008)

-   ANOVA – to evaluate the relationship between recent growth and
    competition for each species and life stage across all sampling
    ewlevations

    -   Had to take the natural log of the response to achieve normality

    -   Used preliminary model selection to compare continuous and
        categorical explan variables

    -   To compare the relative strength of competitive effects on
        adults vs saplings they had to transform the continuous adult
        effects to the categorical index for saplings

<!-- -->

-   Conclusions

    -   Highlights the complexity of range shift dynamics for long-lived
        species ike trees

    -   These results suggest that range contractions will occur due to
        a lag between warming and range shifts

Hillyer & Silman 2010

-   Used a Locally weighted scatterplot smoothing (LOESS) (non-para
    regression technique) to estimate expected survival rates across
    elevation for groups of speices

    -   LOESS doesn’t require a fit function, so is useful for modelling
        processes for which there are no theoretical explanations

Wittich *et al.* 2012

-   Relationship between incident light and Ps described using a
    non-rectangular hyperbolic function

-   ANOVA and scheffes tests compared the means of the three stands

-   Wilcoxon U test used instead of scheffes if data non-normally
    distributed

-   

Ideas for mine

-   Mixed effects models

    -   Random effects

        -   Individual seedling/tree

        -   Time of measurement – to account for differences throughout
            the day in terms of photosynthetic capacity etc.

    -   Fixed effects

        -   


