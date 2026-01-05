import matplotlib.pyplot as plt
import numpy as np

categories = ['Unstratified\nBART', 'Unweighted Logistic\nRegression', 'Stratified\nBART', 'Weighted\nLogistic Regression']

# hard code the values to plot
# values for the total effect
means = [0.23, 0.19, 0.16, 0.18]
lower_ci = [-0.05, 0.05, 0.01, -0.02]
upper_ci = [0.55, 0.36, 0.33, 0.39]

error_lower = np.array(means) - np.array(lower_ci)
error_upper = np.array(upper_ci) - np.array(means)

# values for immediate kidney injury
means_rki = [0.28, 0.04, 0.05, 0.07]
lower_ci_rki = [-0.4, -0.05, -0.29, -0.02]
upper_ci_rki = [0.96, 0.16, 0.4, 0.19]

error_lower_rki = np.array(means_rki) - np.array(lower_ci_rki)
error_upper_rki = np.array(upper_ci_rki) - np.array(means_rki)

# values for all other biological mechanisms
means_obm = [-0.05, 0.15, 0.12, 0.11]
lower_ci_obm = [-0.7, -0.03, -0.23, -0.06]
upper_ci_obm = [0.56, 0.35, 0.41, 0.31]

error_lower_obm = np.array(means_obm) - np.array(lower_ci_obm)
error_upper_obm = np.array(upper_ci_obm) - np.array(means_obm)

x_pos = np.array([1,2,3,4])
x_pos_rki = x_pos + 0.2
x_pos_obm = x_pos + 0.4

plt.figure(figsize=(5,3))

plt.rcParams['font.size'] = 8

plt.errorbar(x_pos, means, yerr=[error_lower, error_upper], fmt='o', capsize=5, color='blue', label='Total Effect')

plt.text(3+0.045, upper_ci[2]+0.05, '$*$', horizontalalignment='right', verticalalignment='center', fontsize=12, color='blue')
plt.text(2+0.045, upper_ci[2]+0.05, '$*$', horizontalalignment='right', verticalalignment='center', fontsize=12, color='blue')
plt.figtext(0.5, 0.01, '* indicates intervals do not cover 0.', wrap=True, horizontalalignment='center', fontsize=6)

# add the label for the immediate kidney injury in the legend
plt.errorbar(x_pos_rki, means_rki, yerr=[error_lower_rki, error_upper_rki], fmt='o', capsize=5, color='purple', label='Immediate Kidney Injury')

# add the label for the all other biological mechanisms in the legend
plt.errorbar(x_pos_obm, means_obm, yerr=[error_lower_obm, error_upper_obm], fmt='o', capsize=5, color='green', label='Other Biological Mechanisms')

# label for the y axis
plt.ylabel('Estimated Causal Risk Difference')

plt.xticks(ticks=x_pos_rki, labels=categories)

plt.axhline(y=0, linestyle='dotted', color='black')

plt.legend()

plt.title('Estimated Total and Decomposed Effects of\nNadir DO$_2$ During XCT and X-Clamp Duration on POAKI')
plt.tight_layout()

# plt.show()
plt.savefig('results.png')
