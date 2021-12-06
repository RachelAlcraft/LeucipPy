
'''
This performs some statistical manipulations on distributions given in series
It has no state
'''

from scipy import stats

def compareDistributionsMannWhitney(distA, distB, sig_level, htmlise=False):
    distA.sort()
    distB.sort()
    u_statistic, p_value = stats.mannwhitneyu(distA,distB)

    diver = "\n"
    if htmlise:
        diver = "<br/>\n"

    test_string = "-- Mann-Whitney U Test --" + diver
    test_string +=  "Hypothesis:The distributions are identical" + diver
    test_string += "Method: If the p-value < " + str(round(sig_level,2)) + " we will reject" + diver
    test_string += "Evidence: The p-value is " + str(round(p_value,4)) + diver
    test_string += "(U Statistic=" + str(round(u_statistic,2)) + ")" + diver
    if p_value < sig_level:
        test_string += "Conclusion: We reject the hypothesis, they are not the same."
    else:
        test_string += "Conclusion: We have no evidence to reject the hypothesis."

    return p_value, test_string


def normalityShapiroWilk(distA, sig_level, htmlise=False):
    distA.sort()
    statistic,pvalue = stats.shapiro(distA)

    diver = "\n"
    if htmlise:
        diver = "<br/>\n"

    test_string = "-- Shapiro-Wilk Normality Test --" + diver
    test_string +=  "Hypothesis:The distribution is normal" + diver
    test_string += "Method: If the p-value < " + str(round(sig_level,2)) + " we will reject" + diver
    test_string += "Evidence: The p-value is " + str(round(pvalue,4)) + diver
    test_string += "(Statistic=" + str(round(statistic,2)) + ")" + diver
    if pvalue < sig_level:
        test_string += "Conclusion: Not normal, we reject the hypothesis."
    else:
        test_string += "Conclusion: Normal, we have no evidence to reject the hypothesis."

    return pvalue, test_string


