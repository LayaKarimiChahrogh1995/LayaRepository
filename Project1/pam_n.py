import numpy as np
from numpy.linalg import matrix_power

text_file = open ('PAM1.txt', 'r')
lines = text_file.readlines()[2:]
n_size = int(input('Enter the unit divergence for PAM matrix \n' ))
def pam_n(n):
    pam_1_matrix = []
    for line in lines:
        pam_1 = line.rstrip().split(',')
        pam_1 = [int(i) for i in pam_1 ] #to int
        pam_1_matrix.append(pam_1)

    pam_arr = np.array(pam_1_matrix)
    #pam_arr = pam_arr.transpose()
    pam_arr = np.true_divide(pam_arr, 10000)
    pam_n_matrix = matrix_power(pam_arr, n)
    print('Probablity mutation matrix is: \n', pam_n_matrix)
    freq = [0.08768333990119, 0.0405129182960021, 0.0408784631518651, 0.0477160345974603, 0.0324709539656211, 0.0378461268859548, 0.0504933695605074, 0.0898249006830963, 0.0328588505954496, 0.0357514442352249, 0.0852464099207531, 0.0791031344407513, 0.0148824394639692, 0.0410010190895639, 0.0515802694709073, 0.0697549720598532, 0.0583275704247605, 0.00931264523877659, 0.0317154088087089, 0.0630397292098708  ]
    #for n in range (1,2):
    #    pam_n = np.multiply(pam_n, arr)
#  finding out score matrix
    sub_matrix_values = np.zeros((20,20))
    for i in range (20):
        for j in range (20):
            sub_matrix_values [i][j] = int(0.5 + 10 * (np.log(pam_n_matrix[j][i] / freq[j]))/np.log(10))
    #print("Substitution score matrix\n", sub_matrix_values)
    return sub_matrix_values

prob_matrix = pam_n(n_size)
print('subs matrix of PAM', n_size , 'is \n', prob_matrix)



