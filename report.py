#
# Plot the report files from the test/test.f90 code.
#
# reads two files: one with all the evaluations (report_all.csv)
# and one with only the best solution found at each evaluation (report_best.csv).
#

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# PREFIX = 'python_'  # process python test files
PREFIX = ''  # process fortran test files.

# file with only the best solution found so far at each evaluation, with columns x1,x2,f
df = pd.read_csv(f'{PREFIX}report_best.csv', header='infer')
x_values = df.iloc[:, 0].values
y_values = df.iloc[:, 1].values
f_values = df.iloc[:, 2].values
i_values = np.arange(len(f_values))

# file with ALL the solutions evaluated, with columns x1,x2,f
df_all = pd.read_csv(f'{PREFIX}report_all.csv', header='infer')
x_values_all = df_all.iloc[:, 0].values
y_values_all = df_all.iloc[:, 1].values
f_values_all = df_all.iloc[:, 2].values
i_values_all = np.arange(len(f_values_all))

##################################################################

# 2D Gradient Plot: All function evaluations with best curve overlay
plt.figure(figsize=(8, 4))

# Create scatter plot of all evaluations with color representing function value
scatter = plt.scatter(x_values_all, y_values_all, c=f_values_all,
                     cmap='coolwarm', alpha=0.6, s=20, edgecolors='none')
plt.colorbar(scatter, label='Function Value $f(x_1, x_2)$')
plt.plot(x_values_all, y_values_all, '-', color='lightgray', linewidth=0.5, alpha=0.6, label='All Evaluations Path', zorder=-1)

# Overlay the best curve (path of best solutions)
plt.plot(x_values, y_values, 'k-', linewidth=2, label='Best Solution Path', zorder=10)
plt.plot(x_values, y_values, 'ko', markersize=4, zorder=11)  # Add markers

# Mark the final best solution
plt.plot(x_values[-1], y_values[-1], 'r*', markersize=10,
         label=f'Final Best: ({x_values[-1]:.3f}, {y_values[-1]:.3f})', zorder=12)

plt.xlabel('$x_1$')
plt.ylabel('$x_2$')
plt.title('2D Function Evaluation Space\n(All Evaluations with Best Path)')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f'{PREFIX}report_2d.jpg', dpi=100)
plt.show()
##################################################################
