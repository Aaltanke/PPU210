# Description:
# Function that performs rainflow counting using 3-point algorithm according to ASTM E1049-85
# Input:
# stress_signal: Row or column vector of stress values from a stress history signal expressed as ndarray.
# Output:
# cycle_counts: Matrix where first column contains cycle counts, second column is the corresponding stress ranges
# and third column contains mean stress values for each range.
# Written by: Louise Åkesson
# Contact: louakes@chalmers.se

# ==========Packages========
import numpy as np
# ==========================


def rainflow_counting(stress_signal):
    """Performs rainflow counting on a given stress history signal with 3-point algorithm.

        :param stress_signal: Vector of stress values from a stress history signal.
        :type stress_signal: numpy.ndarray

        :return: Matrix with cycle counts in first column, stress range in second and mean stress in third.
        :rtype: numpy.ndarray
        """

    # Convert stress_signal to column vector if row vector
    if stress_signal.shape[0] == 1:
        stress_signal = stress_signal.reshape(-1, 1).flatten()

    # Find reversals by checking neighbours of each stress value
    reversals = [stress_signal[0]]

    for stress in range(1, len(stress_signal) - 1):
        if stress_signal[stress] > stress_signal[stress - 1] and \
                stress_signal[stress] > stress_signal[stress + 1]:  # Local maximum (peak)
            reversals.append((stress_signal[stress]))
        elif stress_signal[stress] < stress_signal[stress - 1] and \
                stress_signal[stress] < stress_signal[stress + 1]:  # Local minimum (valley)
            reversals.append((stress_signal[stress]))

    reversals.append(stress_signal[-1])

    # Initialize cycle counts and current points
    counts = []
    current_points = []

    # Step 1: Read the next reversal. If out of data, go to Step 6.
    for reversal in range(len(reversals)):

        current_points.append(reversal)
        number_of_points = len(current_points)

        # Step 2: Less than 3 points - go to step 1. Define X and Y from 3 latest reversals that are not discarded
        while number_of_points >= 3:
            x_range = abs(reversals[current_points[number_of_points - 2]] -
                          reversals[current_points[number_of_points - 1]])  # Calculate X range
            y_range = abs(reversals[current_points[number_of_points - 3]] -
                          reversals[current_points[number_of_points - 2]])  # Calculate Y range

            # Step 3: Compare absolute values of X and Y
            if x_range >= y_range:  # Step 3b): If X>=Y, go to Step 4.
                y_mean = 0.5 * (reversals[current_points[number_of_points - 3]] +
                                reversals[current_points[number_of_points - 2]])  # Calculate cycle mean

                # Step 4: If Y contains starting point S (always when current points are 3), go to Step 5
                if number_of_points == 3:

                    # Step 5: Count as half cycle and discard first point in Y. Go to Step 2
                    counts.append([0.5, y_range, y_mean])
                    del current_points[0]
                    number_of_points = len(current_points)

                # Step 4: If Y does not contain S, count as one cycle, discard both points of Y and go to Step 2
                else:
                    counts.append([1.0, y_range, y_mean])
                    del current_points[number_of_points - 3]
                    del current_points[number_of_points - 3]
                    number_of_points = len(current_points)

            # Step 3b): If X<Y, go to Step 1
            else:
                break

    # Step 6: Count remaining reversals as half cycles
    for residual in range(number_of_points - 1):
        y_range = abs(reversals[current_points[residual]] - reversals[current_points[residual + 1]])
        y_mean = 0.5 * (reversals[current_points[residual]] + reversals[current_points[residual + 1]])
        counts.append([0.5, y_range, y_mean])

    cycle_counts = np.array(counts)
    
    # Histogram parameters
    hist, bin_edges = np.histogram(cycle_counts[:, 1], bins=65, weights=cycle_counts[:, 0])

    return cycle_counts, hist, bin_edges     
  