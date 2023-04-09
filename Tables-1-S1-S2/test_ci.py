from p_from_scaled_containment import *
from mutation_model_simulator import MutationModel
import numpy as np
import sys

def progressbar(it, prefix="", size=60, file=sys.stdout):
    count = len(it)
    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()

if __name__ == '__main__':
    mutation_rates = [0.001, 0.1, 0.2]
    scale_factors = [0.1, 0.2, 0.05]
    confidence = 0.95
    k_mer_lengths = [21, 31, 51, 71, 100]
    num_k_mers_list = [10000, 100000, 1000000]
    num_simulations = 10000
    '''
    # this code segment generates variance of C_scale using multiple simulations
    for k_mer_length in k_mer_lengths:
        mutation_model = MutationModel(num_k_mers+k_mer_length-1, k_mer_length, mutation_rate, scale_factor)
        scaled_containment_indices = []
        for iter_count in progressbar(range(num_simulations), " Simulating: ", 40):
            mutation_model.generate()
            scaled_containment_indices.append(mutation_model.count_scaled_containment())
        print(np.var(scaled_containment_indices))
    '''
    #ci_list = compute_ci_from_jaccard_test([0.248439450687], num_k_mers, 21, confidence, scale_factor)
    #print(ci_list)

    #print('OK')

    #this segment generates scaled jaccard and keeps it in a list
    '''
    for k_mer_length in k_mer_lengths:
        mutation_model = MutationModel(num_k_mers+k_mer_length-1, k_mer_length, mutation_rate, scale_factor)
        scaled_jaccard_indices = []
        for iter_count in range(num_simulations):
            mutation_model.generate()
            scaled_jaccard_indices.append(mutation_model.count_scaled_jaccard())
        ci_from_scaled_jaccard = compute_ci_from_jaccard(scaled_jaccard_indices, num_k_mers,
                                                         k_mer_length, confidence, scale_factor)
        hit_count = 0
        miss_count = 0
        for ci in ci_from_scaled_jaccard:
            plow = ci[4]
            phigh = ci[5]
            if mutation_rate >= plow and mutation_rate <= phigh:
                hit_count += 1
            else:
                miss_count += 1
            #print(mutation_rate, plow, phigh, ci[3])
        print(k_mer_length, 1.0*hit_count/(hit_count+miss_count), hit_count, miss_count)

        ci_from_scaled_jaccard = compute_ci_from_jaccard_test(scaled_jaccard_indices, num_k_mers,
                                                         k_mer_length, confidence, scale_factor)
        hit_count = 0
        miss_count = 0
        for ci in ci_from_scaled_jaccard:
            plow = ci[4]
            phigh = ci[5]
            if mutation_rate >= plow and mutation_rate <= phigh:
                hit_count += 1
            else:
                miss_count += 1
            #print(mutation_rate, plow, phigh, ci[3])
        print(k_mer_length, 1.0*hit_count/(hit_count+miss_count), hit_count, miss_count)
    '''

    '''
    # this code segment generates % of times mut_rate is with in the conf.interval that we calculate
    for k_mer_length in k_mer_lengths:
        mutation_model = MutationModel(num_k_mers+k_mer_length-1, k_mer_length, mutation_rate, scale_factor)
        scaled_containment_indices = []
        for iter_count in range(num_simulations):
            mutation_model.generate()
            scaled_containment_indices.append(mutation_model.count_scaled_containment())
        ci_list = compute_confidence_interval_test3(scaled_containment_indices, num_k_mers,
                                                                            k_mer_length,
                                                                            confidence,
                                                                            scale_factor)
        hit_count = 0
        miss_count = 0
        for ci in ci_list:
            plow = ci[4]
            phigh = ci[5]
            if mutation_rate >= plow and mutation_rate <= phigh:
                hit_count += 1
            else:
                miss_count += 1
            #print(mutation_rate, plow, phigh)
        print(k_mer_length, 1.0*hit_count/(hit_count+miss_count), hit_count, miss_count)

        # this code segment will generate CI in one step and check what % the true mut_rate falls with in our conf.interval
        #ci_one_step_list = compute_confidence_interval_one_step(scaled_containment_indices, num_k_mers,
        #                                                        k_mer_length, confidence, scale_factor)

        ci_one_step_list = compute_confidence_interval_test(scaled_containment_indices, num_k_mers,
                                                                k_mer_length, confidence, scale_factor)

        hit_count = 0
        miss_count = 0
        for ci in ci_one_step_list:
            C_scale = ci[3]
            plow = ci[4]
            phigh = ci[5]
            #print(C_scale, mutation_rate, plow, phigh)
            if mutation_rate >= plow and mutation_rate <= phigh:
                hit_count += 1
            else:
                miss_count += 1
        print(k_mer_length, 1.0*hit_count/(hit_count+miss_count), hit_count, miss_count)
        '''
    for scale_factor in scale_factors:
        max_trials = int(1.0/scale_factor)+1
        for num_k_mers in num_k_mers_list:
            for mutation_rate in mutation_rates:
                for k_mer_length in k_mer_lengths:
                    mutation_model = MutationModel(num_k_mers+k_mer_length-1, k_mer_length, mutation_rate, scale_factor)
                    scaled_containment_indices = []
                    for iter_count in range(num_simulations):
                        counter = 0
                        while counter < max_trials:
                            mutation_model.generate()
                            candidate_scaled_contaiment_index = mutation_model.count_scaled_containment()
                            if candidate_scaled_contaiment_index == 0.0 or candidate_scaled_contaiment_index == 1.0:
                                counter += 1
                            else:
                                break
                        scaled_containment_indices.append(candidate_scaled_contaiment_index)
                    ci_one_step_list = compute_confidence_interval_one_step(scaled_containment_indices, num_k_mers,
                                                                            k_mer_length, confidence, scale_factor)

                    hit_count = 0
                    miss_count = 0
                    above_mid_point_count = 0
                    below_mid_point_count = 0
                    for ci in ci_one_step_list:
                        C_scale = ci[3]
                        plow = ci[4]
                        phigh = ci[5]
                        p_point = ci[6]
                        #print(C_scale, mutation_rate, plow, phigh)
                        if mutation_rate >= plow and mutation_rate <= phigh:
                            hit_count += 1
                        else:
                            miss_count += 1
                        mid_point = (plow + phigh) / 2.0
                        if mutation_rate >= mid_point:
                            above_mid_point_count += 1
                        else:
                            below_mid_point_count += 1
                    print(k_mer_length, mutation_rate, num_k_mers, scale_factor, 1.0*hit_count/(hit_count+miss_count), hit_count, miss_count, above_mid_point_count, below_mid_point_count)
    '''
                    ci_two_step_list = compute_confidence_interval_test3(scaled_containment_indices, num_k_mers,
                                                                            k_mer_length, confidence, scale_factor)
                    hit_count = 0
                    miss_count = 0
                    above_mid_point_count = 0
                    below_mid_point_count = 0
                    for ci in ci_two_step_list:
                        C_scale = ci[3]
                        plow = ci[4]
                        phigh = ci[5]
                        #p_point = ci[6]
                        #print(C_scale, mutation_rate, plow, phigh)
                        if mutation_rate >= plow and mutation_rate <= phigh:
                            hit_count += 1
                        else:
                            miss_count += 1
                        mid_point = (plow + phigh) / 2.0
                        if mutation_rate >= mid_point:
                            above_mid_point_count += 1
                        else:
                            below_mid_point_count += 1
                    print(k_mer_length, mutation_rate, num_k_mers, 1.0*hit_count/(hit_count+miss_count), hit_count, miss_count, above_mid_point_count, below_mid_point_count)
            #'''
