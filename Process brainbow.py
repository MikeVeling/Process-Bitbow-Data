#################################################################################
#                            User editable varables                             #
#################################################################################
include_larva_data_conditions=False                                             # Boolean allowing for a by_the_larva output collecting all larval level data in the input file level output
include_larva_data_dates=False                                                  # Boolean allowing for a by_the_larva output collecting all larval level data in the date level output
include_by_the_larva_outputs=False                                              # Boolean allowing for a by_the_larva outputs for the color_stat, group_color_stat, and color_number_stat output files
include_pvalue_conditions=True                                                  # Boolean allowing for the p-value calculation at the input file level output
include_pvalue_dates=False                                                      # Boolean allowing for the p-value calculation at the date level output (warning, computationally intensive)
include_pvalue_larva=False                                                      # Boolean allowing for the p-value calculation at the larva level output (warning, computationally intensive)
include_conditions_in_flaten=False                                              # Boolean allowing for the generation of a set of flat files containing all larva level data in the input file level output
include_date_in_flaten=False                                                    # Boolean allowing for the generation of a set of flat files containing all larva level data in the date level output
include_larva_in_flaten=False                                                   # Boolean allowing for the generation of a set of flat files containing all larva level data in the larva level output
include_randomized_match_data=True                                              # Boolean allowing for the inclusion of randomized match data in the by_the_colors.csv output. This is useful for comparing match numbers to the number of times neurons randomly matched
random_only_colors=True                                                         # Boolean determining if only neurons with color are randomized
test_stat='sort_of_probibility_of_matches'                                      # This can be 'number_of_matches' or 'sort_of_probibility_of_matches' or 'sort_of_prob_of_matches_times_misses' or 'percent_match'
Multiple_hypothesis_method='fdr_bh'                                             # This var inherits its testing method from statsmodels.stats.multitest.multipletests
alpha_level=0.05                                                                # This var sets the alpha level for statsmodels.stats.multitest.multipletests
number_of_p_val_itterations=100                                                 # This is the number of iterations used to generate the background distribution. It is set to 100 so the code can run on most computers, but I recommend at least 100,000 for real analysis. 
itteration_print=100                                                            # How many iterations do you want to run per thread. With my i7 5820 12 core 100 was about the right number for about 4 to 6 min of processing time
Thread_multiplyer=1                                                             # how many simulations jobs do you want to run based on the number of threads your computer has. Default is 1 meaning it will use all of your threads.
Thread_offset=0                                                                 # how many extra jobs do you want to run above the number of cores you have. This can also be negative to allow for running of other programs while this is running.
multithread_sleep_time=1                                                        # how often should python check for finished jobs
verbose_multithread_worker_search=True                                          # Should python tell you while it is searching for idle workers
only_include_some_segments=True                                                 # When False, this will simply perform all the hemisegment analysis
path_of_key_file='key.csv'                                                      # Pathway of the key.csv file
included_hemisegments_path='Active_hemisegments.csv'                            # pathway of the active_hemisegments.csv file
exclude_color_layers=['Red','Yellow','Blue','Green','Purple']                   # These are neuron grouping layers that should not be used in the analysis. This allows you to keep a complex key file but maintain a simple analysis
#################################################################################
#                             Defineing my objects                              #
#################################################################################
import warnings, itertools, sys, json                                           #
from statsmodels.stats.multitest import multipletests                           #
from operator import itemgetter                                                 #
from numpy import std as stdev                                                  #
class basefunctions:                                                            #
    def get_larva_list(self):                                                   #
        if isinstance(self,larva):                                              #
            larva_objects_list=[self]                                           #
        elif isinstance(self,date):                                             #
            larva_objects_list=self.larvae                                      #
        elif isinstance(self,condition):                                        #
            larva_objects_list=[]                                               #
            for date_object in self.dates:                                      #
                for larva_object in date_object.larvae:                         #
                    larva_objects_list.append(larva_object)                     #
        return larva_objects_list                                               #
    def get_hemisegments_list(self):                                            #
        hemisegments_list=[]                                                    #
        for larva_object in self.get_larva_list():                              #
            for hemisegment_object in  larva_object.hemisegments:               #
                hemisegments_list.append(hemisegment_object)                    #  
        return hemisegments_list                                                #
    def get_flat_hemisegment_table(self):                                       #
        flat_data_table=[]                                                      #
        for hemisegment_object in self.get_hemisegments_list():                 #
            flat_data_table.append(hemisegment_object.neuron_data)              #
        self.flat_data_table=flat_data_table                                    #
        return flat_data_table                                                  #
    def get_flat_color_list(self):                                              # Here, I define functions to get flat data about the colors observed. It also updates universal numbers used in the dataset like number of neurons observed and number of neurons with color
        self.number_of_neurons_not_observed=0                                   #
        self.total_neurons_observed=0                                           #
        self.total_neurons_with_any_color=0                                     #
        flat_data=[]                                                            #
        data_table=self.get_flat_hemisegment_table()                            #
        for row in data_table:                                                  #
            for col in row:                                                     #
                flat_data.append(col)                                           #
        for row in flat_data:                                                   #
            if row == '':                                                       #
                self.number_of_neurons_not_observed+=1                          #
            else:                                                               #
                self.total_neurons_observed+=1                                  #
                active_chanels=sum_active_chanels(row)                          #
                if active_chanels != 0:                                         #
                    self.total_neurons_with_any_color+=1                        #
        return flat_data                                                        #
    def make_color_number_stats(self):                                          # This function uses the flat data generated ablve and makes the color_number_stats output
        flat_data=self.get_flat_color_list()                                    #
        output_file=[['total number of chanels with color',                     #
                      'Total with given number of colors',                      #
                      'Total neurons',                                          #
                      'Total neurons with any color']]                          #
        neuron_number_count_dic={}                                              #
        for color_count in range(0,6):                                          #
            neuron_number_count_dic[color_count]=0                              #
        for row in flat_data:                                                   #
            if row != '':                                                       #
                active_chanels=sum_active_chanels(row)                          #
                neuron_number_count_dic[active_chanels]+=1                      #
        for active_chanel_count in neuron_number_count_dic:                     #
            next_line=[str(active_chanel_count),                                #
                       str(neuron_number_count_dic[active_chanel_count]),       #
                       str(self.total_neurons_observed),                        #
                       str(self.total_neurons_with_any_color)]                  #
            output_file.append(next_line)                                       #
        return output_file                                                      #
    def make_color_stats(self):                                                 # This function uses the flat data generated ablve and makes the color_stats output
        flat_data=self.get_flat_color_list()                                    #
        output_file=[['',                                                       #
                      'Color chanel 1',                                         #
                      'Color chanel 2',                                         #
                      'Color chanel 3',                                         #
                      'Color chanel 4',                                         #
                      'Color chanel 5']]                                        #
        color_activ_dic={}                                                      #
        for chanel_ID in range(0,5):                                            #
            color_activ_dic[chanel_ID]=0                                        #
        for datium in flat_data:                                                #
            if datium != '':                                                    #
                for datium_col in color_activ_dic:                              #
                    point=datium[datium_col]                                    #
                    assert point == '1' or point == '0'                         #
                    if point == '1':                                            #
                        color_activ_dic[datium_col]+=1                          #
        next_line=['count with chanel',                                         #
                    str(color_activ_dic[0]),                                    #
                    str(color_activ_dic[1]),                                    #
                    str(color_activ_dic[2]),                                    #
                    str(color_activ_dic[3]),                                    #
                    str(color_activ_dic[4])]                                    #
        output_file.append(next_line)                                           #
        next_line=['Total Neurons',                                             #
                   str(self.total_neurons_observed),                            #
                   str(self.total_neurons_observed),                            #
                   str(self.total_neurons_observed),                            #
                   str(self.total_neurons_observed),                            #
                   str(self.total_neurons_observed)]                            #
        output_file.append(next_line)                                           #
        next_line=['Total Neurons with color',                                  #
                   str(self.total_neurons_with_any_color),                      #
                   str(self.total_neurons_with_any_color),                      #
                   str(self.total_neurons_with_any_color),                      #
                   str(self.total_neurons_with_any_color),                      #
                   str(self.total_neurons_with_any_color)]                      #
        output_file.append(next_line)                                           #
        return output_file                                                      #
    def make_group_color_stats(self):                                           # This function uses the flat data generated ablve and makes the color_group_stats output
        if hasattr(self, 'group_color_stat_output'):                            #
            return self.group_color_stat_output                                 #
        flat_data=self.get_flat_color_list()                                    #
        colors_list=make_all_colors_order()                                     #
        self.color_match_prob_dic={}                                            #
        self.color_match_prob_list=[]                                           #
        output=[['Color ID',                                                    #
                 'Total number of chanels with color',                          #
                 'Total with given color',                                      #
                 'Total neurons',                                               #
                 'Total neurons with any color',                                #
                 'Probability of the given color (including blanks)',           #
                 'Probability of the given color (excluding blanks)',           #
                 'Total posible colors',                                        #
                 'Total posible non-blank colors']]                             #
        colors_dic={}                                                           #
        for color in colors_list:                                               #
            colors_dic[color]=0                                                 #
        for color in flat_data:                                                 #
            if color != '':                                                     #
                colors_dic[color]+=1                                            #
        self.prob_of_match_including_blanks=0.0                                 #
        self.prob_of_match_excluding_blanks=0.0                                 #
        for row in colors_list:                                                 #
            count_of_active_chanels=str(sum_active_chanels(row))                #
            if self.total_neurons_with_any_color != 0 and row!='00000':         #
                excluding_blanks_prob=(float(colors_dic[row])/                  #
                                       float(self.total_neurons_with_any_color))#
                self.prob_of_match_excluding_blanks+=excluding_blanks_prob**2   #
                if excluding_blanks_prob != 0.0:                                #
                    self.color_match_prob_dic[row]=excluding_blanks_prob**2     #
                    self.color_match_prob_list.append([row,                     #
                                                       excluding_blanks_prob])  #
            else:                                                               #
                excluding_blanks_prob='N/A'                                     #
            including_blanks_prob=(float(colors_dic[row])/                      #
                                   float(self.total_neurons_observed))          #
            if row != '00000':                                                  #
                self.prob_of_match_including_blanks+=including_blanks_prob**2   #
            next_line=['"'+row+'"',                                             #
                       count_of_active_chanels,                                 #
                       str(colors_dic[row]),                                    #
                       str(self.total_neurons_observed),                        #
                       str(self.total_neurons_with_any_color),                  #
                       str(including_blanks_prob),                              #
                       str(excluding_blanks_prob),                              #
                       str(32),                                                 #
                       str(31)]                                                 #
            output.append(next_line)                                            #
        output[0]+=['',                                                         #
                    'Probability of a match (including blanks)',                #
                    'Probability of a match (excluding blanks)']                #
        output[1]+=['',                                                         #
                    str(self.prob_of_match_including_blanks),                   #
                    str(self.prob_of_match_excluding_blanks)]                   #
        self.color_match_prob_list=sorted(self.color_match_prob_list,           #
                                          key=itemgetter(1))[::-1]              #
        self.simp_color_match_prob_list=[]                                      #
        for row in self.color_match_prob_list:                                  #
            self.simp_color_match_prob_list.append(row[0])                      #
        self.group_color_stat_output=output                                     #
        return output                                                           #
    def make_match_report(self,                                                 #
                          neuron_1,                                             #
                          neuron_2,                                             #
                          include_header=False,                                 #
                          include_pvalues=False,                                #
                          pval_dic_path='N/A'):                                 #
        assert isinstance(self.prob_of_match_including_blanks, float)           #
        assert isinstance(self.prob_of_match_excluding_blanks, float)           #
        assert isinstance(neuron_1, neuron)                                     #
        assert isinstance(neuron_2, neuron)                                     #
        assert neuron_1.layer == neuron_2.layer                                 #
        neuron_layer=neuron_1.layer                                             #
        if neuron_1 == neuron_2:                                                #
            rep_neurons=True                                                    #
        else:                                                                   #
            rep_neurons=False                                                   #
        output_color_order=make_all_colors_order()                              #
        output_color_order_header=[]                                            #
        for row in output_color_order:                                          #
            output_color_order_header.append('"'+row+'"')                       #
        header=['Neuron 1',                                                     #
                'Neuron 2',                                                     #
                'Neuron 1 group size',                                          #
                'Neuron 2 group size',                                          #
                'Simple Match Name',                                            #
                'Total number of times both neurons were observed',             #
                'Total number of times both neurons were observed '+            #
                'and at least one had a color',                                 #
                'Total number of times both neurons had color',                 #
                'Total number of times both neurons had the same color '+       #
                '(matches)',                                                    #
                'Total number of times both neurons were observed and at '+     #
                'least one had a color - Total number of times both neurons '+  #
                'had the same color (non-matches)',                             #
                'Fractional match']                                             #
        if include_pvalues and neuron_layer not in exclude_color_layers:        #
            try:                                                                #
                working_pval_dic=self.pval_dic[neuron_layer]                    #
            except:                                                             #
                working_pval_dic=self.make_pval_dic(number_of_p_val_itterations,#
                                                   neuron_layer,                #
                                                   neurons_dic,                 #
                                                   itteration_print,            #
                                                   pval_dic_path)               #
                self.pval_dic[neuron_layer]=working_pval_dic                    #
            header+=['Relatedness test statistic']                              #
            header+=['Relatedness test P-Value']                                #
            if include_randomized_match_data:                                   #
                header+=['Total number of times both neurons were observed '+   #
                        'and at least one had a color randomly']                #
                header+=['Total number of times both neurons had the same '+    #
                         'color randomly (matches)']                            #
                header+=['Total number of times both neurons were observed '+   #
                         'and at least one had a color randomly - '+            #
                         'Total number of times both neurons had the same '+    #
                         'color randomly (non-matches)']                        #
        header+=['']+output_color_order_header                                  #
        neuron_1_string=neuron_1.identifier                                     #
        neuron_2_string=neuron_2.identifier                                     #
        neuron_1_group_size=len(neuron_1.col_IDs)                               #
        neuron_2_group_size=len(neuron_2.col_IDs)                               #
        simple_match_name=neuron_1_string+' to '+neuron_2_string                #
        matches_dic={}                                                          #
        for color in output_color_order:                                        #
            matches_dic[color]=0                                                #
        both_neruons_observed=0                                                 #
        both_neurons_observed_one_with_color=0                                  #
        both_neurons_with_color=0                                               #
        matches=0                                                               #
        for hemisegment_obj in self.get_hemisegments_list():                    #
            working_neuron_data=hemisegment_obj.neuron_data                     #
            neuron_1_colors=[]                                                  #
            neuron_1_coordanates=neuron_1.col_IDs                               #
            for neuron_1_coordanate in neuron_1_coordanates:                    #
                neuron_1_colors.append(working_neuron_data[neuron_1_coordanate])#
            neuron_2_colors=[]                                                  #
            neuron_2_coordanates=neuron_2.col_IDs                               #
            for neuron_2_coordanate in neuron_2_coordanates:                    #
                neuron_2_colors.append(working_neuron_data[neuron_2_coordanate])#
            if both_have_data(neuron_1_colors,                                  #
                              neuron_2_colors,                                  #
                              rep_neurons=rep_neurons):                         #
                both_neruons_observed+=1                                        #
                if either_have_color(neuron_1_colors,                           #
                                     neuron_2_colors,                           #
                                     rep_neurons=rep_neurons):                  #
                    both_neurons_observed_one_with_color+=1                     #
                    if both_have_color(neuron_1_colors,                         #
                                       neuron_2_colors,                         #
                                       rep_neurons=rep_neurons):                #
                        both_neurons_with_color+=1                              #
                        if both_have_match(neuron_1_colors,                     #
                                           neuron_2_colors,                     #
                                           rep_neurons=rep_neurons):            #
                            matches+=1                                          #
                            matching_colors_list=find_matches(neuron_1_colors,  #
                                                              neuron_2_colors,  #
                                                        rep_neurons=rep_neurons)#
                            for color_match in matching_colors_list:            #
                                matches_dic[color_match]+=1                     #
        color_matches_list=[]                                                   #
        for output_color in output_color_order:                                 #
            color_matches_list.append(matches_dic[output_color])                #
        if float(both_neurons_observed_one_with_color) != float(0):             #
            percent_match=str(float(matches)/                                   #
                              float(both_neurons_observed_one_with_color))      #
        else:                                                                   #
            percent_match='N/A'                                                 #
        non_matches=both_neurons_observed_one_with_color-matches                #
        next_line=[neuron_1_string,                                             #
                   neuron_2_string,                                             #
                   neuron_1_group_size,                                         #
                   neuron_2_group_size,                                         #
                   simple_match_name,                                           #
                   both_neruons_observed,                                       #
                   both_neurons_observed_one_with_color,                        #
                   both_neurons_with_color,                                     #
                   matches,                                                     #
                   non_matches,                                                 #
                   percent_match]                                               #
        if include_pvalues:                                                     #
            if simple_match_name in working_pval_dic:                           #
                proc_simple_match_name=simple_match_name                        #
            else:                                                               #
                simple_match_name_rev_1=simple_match_name.split(' to ')[0]      #
                simple_match_name_rev_2=simple_match_name.split(' to ')[1]      #
                proc_simple_match_name=(simple_match_name_rev_2+' to '+         #
                                       simple_match_name_rev_1)                 #
            pval_conversion_dic=working_pval_dic[proc_simple_match_name]        #
            random_total=pval_conversion_dic['totals']                          #
            random_matches=pval_conversion_dic['matches']                       #
            test_stat_val,p_val=self.get_test_stat_and_pval(neuron_1,           #
                                                            neuron_2,           #
                                                            pval_conversion_dic,#
                                                            neurons_dic,        #
                                                            neuron_layer)       #
            next_line+=[test_stat_val]                                          #
            next_line+=[p_val]                                                  #
            if include_randomized_match_data:                                   #
                next_line+=[random_total]                                       #
                next_line+=[random_matches]                                     #
                next_line+=[random_total-random_matches]                        #
        next_line+=['']+color_matches_list                                      #
        if include_header:                                                      #
            return[header,next_line]                                            #
        elif not include_header:                                                #
            return next_line                                                    #
    def make_match_report_larva(self,neuron_1,neuron_2):                        #
        larva_list=self.get_larva_list()                                        #
        output_list=[]                                                          #
        for larva_object in larva_list:                                         #
            larva_info=larva_object.larva.split('|||')                          #
            if output_list == []:                                               #
                output_data_pre=larva_object.make_match_report(neuron_1,        #
                                                               neuron_2,        #
                                                            include_header=True)#
                output_list.append(['Condition',                                #
                                    'Date',                                     #
                                    'Larva']+output_data_pre[0])                #
                output_data=output_data_pre[1]                                  #
            else:                                                               #
                output_data=larva_object.make_match_report(neuron_1,neuron_2)   #
            output_list.append(larva_info+output_data)                          #
        return output_list                                                      #
    def get_test_stat_and_pval(self,                                            #
                               neuron_1,                                        #
                               neuron_2,                                        #
                               pval_conversion_dic,                             #
                               neurons_dic,neuron_layer):                       #
        neuron_1_ID=neuron_1.identifier                                         #
        neuron_2_ID=neuron_2.identifier                                         #
        try:                                                                    #
            self.flat_data_table                                                #
        except:                                                                 #
            self.get_flat_hemisegment_table()                                   #
        input_table=self.flat_data_table                                        #
        color_match_prob_dic=self.color_match_prob_dic                          #
        coords_1=neurons_dic[neuron_layer][neuron_1_ID].col_IDs                 #
        coords_2=neurons_dic[neuron_layer][neuron_2_ID].col_IDs                 #
        if include_randomized_match_data:                                       #
            [test_stat_value,                                                   #
            both_neurons_observed_one_with_color,                               #
            matches]=get_test_stat(input_table,                                 #
                                   coords_1,                                    #
                                   coords_2,                                    #
                                   color_match_prob_dic)                        #
        else:                                                                   #
            test_stat_value=get_test_stat(input_table,                          #
                                          coords_1,                             #
                                          coords_2,                             #
                                          color_match_prob_dic)                 #
        if test_stat_value in pval_conversion_dic:                              #
            proc_test_stat_val=test_stat_value                                  #
        else:                                                                   #
            proc_test_stat_val='something_random'                               #
            keys=only_numb_keys(pval_conversion_dic)                            #
            if test_stat_value < min(keys):                                     #
                proc_test_stat_val=min(keys)                                    #
            elif test_stat_value > max(keys):                                   #
                proc_test_stat_val=max(keys)                                    #
            else:                                                               #
                for key_id in range(0,len(keys)):                               #
                    if (test_stat_value > keys[key_id] and                      #
                            test_stat_value < keys[key_id+1]):                  #
                            if (test_stat=='number_of_matches' or               #
                                test_stat=='percent_match'):                    #
                                proc_test_stat_val=keys[key_id]                 #
                            elif (test_stat=='sort_of_probibility_of_matches' or#
                             test_stat=='sort_of_prob_of_matches_times_misses'):#
                                proc_test_stat_val=keys[key_id+1]               #
            if proc_test_stat_val == 'something_random':                        #
                print(test_stat_value)                                          #
                print(keys)                                                     #
                assert proc_test_stat_val != 'something_random'                 #
        p_val=pval_conversion_dic[proc_test_stat_val]                           #
        return test_stat_value, p_val                                           #
    def make_color_number_stats_larva(self):                                    #
        larva_list=self.get_larva_list()                                        #
        return_list=[]                                                          #
        for larva_obj in larva_list:                                            #
            return_list.append(larva_obj.larva.split('|||'))                    #
            for row in larva_obj.make_color_number_stats():                     #
                return_list.append(row)                                         #
        return return_list                                                      #
    def make_color_stats_larva(self):                                           #
        larva_list=self.get_larva_list()                                        #
        return_list=[]                                                          #
        for larva_obj in larva_list:                                            #
            return_list.append(larva_obj.larva.split('|||'))                    #
            for row in larva_obj.make_color_stats():                            #
                return_list.append(row)                                         #
        return return_list                                                      #
    def make_group_color_stats_larva(self):                                     #
        larva_list=self.get_larva_list()                                        #
        return_list=[]                                                          #
        for larva_obj in larva_list:                                            #
            return_list.append(larva_obj.larva.split('|||'))                    #
            for row in larva_obj.make_group_color_stats():                      #
                return_list.append(row)                                         #
        return return_list                                                      #
    def make_random_table(self):                                                #
        try:                                                                    #
            self.flat_data_table                                                #
        except:                                                                 #
            self.get_flat_hemisegment_table()                                   #
        return randomize_objects_in_table(self.flat_data_table)                 #
    def make_pval_dic(self,                                                     #
                      number_of_p_val_itterations,                              #
                      neuron_layer,                                             #
                      neurons_dic,                                              #
                      itteration_print,                                         #
                      pval_dic_path):                                           #
        print('\ton the '+neuron_layer+' layer.')                               #
        output_pval_path=pval_dic_path+'pval_dic.json'                          #
        self.make_group_color_stats()                                           #
        if os.path.isfile(output_pval_path):                                    #
            simp_dic=load_pval_dic(output_pval_path)                            #
            if simp_dic['itterations'] != number_of_p_val_itterations:          #
                os.remove(output_pval_path)                                     #
            else:                                                               #
                return simp_dic                                                 #
        if not number_of_p_val_itterations % itteration_print == 0:             #
            print('your number of itterations is not divisable by your print '+ #
                  'itterations. Fix this')                                      #
        assert number_of_p_val_itterations % itteration_print == 0              #
        number_of_jobs=int(number_of_p_val_itterations/itteration_print)        #
        final_dic=make_test_stat_dic_multithreaded(self,                        #
                                                   itteration_print,            #
                                                   neuron_layer,                #
                                                   neurons_dic,                 #
                                                   neuron_pairs_name_dic,       #
                                                   number_of_jobs,              #
                                                   pval_dic_path)               #
        json_text_dic=json.dumps(final_dic)                                     #
        open(output_pval_path,'w').write(json_text_dic)                         #
        return final_dic                                                        #
class condition(basefunctions):                                                 #
    def __init__(self,fine_name):                                               #
        assert isinstance(fine_name, basestring)                                #
        self.file_name=fine_name                                                #
        self.name=fine_name                                                     #
        self.dates=[]                                                           #
        self.type='contition'                                                   #
        self.pval_dic={}                                                        #
    def add_date_relationship(self,date_obj):                                   #
        assert isinstance(date_obj, date)                                       #
        if date_obj not in self.dates:                                          #
            self.dates.append(date_obj)                                         #
        else:                                                                   #
            warnings.warn('It looks like your date was already in your '        #
                          'condition, just so you know',SyntaxWarning)          #
class date(basefunctions):                                                      #
    def __init__(self,date_str,condition_obj):                                  #
        assert isinstance(date_str, basestring)                                 #
        self.date=date_str                                                      #
        self.name=date_str                                                      #
        assert isinstance(condition_obj,condition)                              #
        self.condition=condition_obj                                            #
        condition_obj.add_date_relationship(self)                               #
        self.larvae=[]                                                          #
        self.type='date'                                                        #
        self.pval_dic={}                                                        #
    def add_larva_relationship(self,larva_obj):                                 #
        assert isinstance(larva_obj, larva)                                     #
        if larva_obj not in self.larvae:                                        #
            self.larvae.append(larva_obj)                                       #
        else:                                                                   #
            warnings.warn('It looks like your larva was already in your '       #
                          'date object, just so you know.',SyntaxWarning)       #
class larva(basefunctions):                                                     #
    def __init__(self,larva_str,date_obj):                                      #
        assert isinstance(larva_str, basestring)                                #
        self.larva=larva_str                                                    #
        self.name=larva_str                                                     #
        assert isinstance(date_obj,date)                                        #
        self.date=date_obj                                                      #
        date_obj.add_larva_relationship(self)                                   #
        self.hemisegments=[]                                                    #
        self.type='larva'                                                       #
        self.pval_dic={}                                                        #
    def add_hemisegment_relationship(self,hemisegment_obj):                     #
        assert isinstance(hemisegment_obj, hemisegment)                         #
        if hemisegment_obj not in self.hemisegments:                            #
            self.hemisegments.append(hemisegment_obj)                           #
        else:                                                                   #
            warnings.warn('It looks like your hemisegment was already in your ' #
                          'larva object, just so you know.',SyntaxWarning)      #
class hemisegment:                                                              #
    def __init__(self,hemisegment_str,larva_obj):                               #
        assert isinstance(hemisegment_str, basestring)                          #
        self.hemisegment=hemisegment_str                                        #
        assert isinstance(larva_obj,larva)                                      #
        self.larva=larva_obj                                                    #
        larva_obj.add_hemisegment_relationship(self)                            #
    def add_data(self,neuron_data):                                             #
        self.neuron_data=neuron_data                                            #
class neuron:                                                                   #
    def __init__(self,neuron_identifier,layer):                                 #
        assert isinstance(hemisegment_str, basestring)                          #
        self.identifier=neuron_identifier                                       #
        self.layer=layer                                                        #
        self.col_IDs=[]                                                         #
    def add_col_ID (self,ID):                                                   #
        assert isinstance(ID, int)                                              #
        self.col_IDs.append(ID)                                                 #
#################################################################################
#                     Defineing my data processing functions                    #
#################################################################################
from random import shuffle                                                      #
import collections                                                              #
def float_list_special(lst):                                                    #
    return_list=[]                                                              #
    totals=0                                                                    #
    matches=0                                                                   #
    for row_pre in lst:                                                         #
        if row_pre != '':                                                       #
            if include_randomized_match_data:                                   #
                [data_pre,totals_trial,matches_trial]=row_pre.split('|||')      #
                totals+=int(totals_trial)                                       #
                matches+=int(matches_trial)                                     #
            else:                                                               #
                data_pre=row_pre                                                #
            if data_pre != 'N/A':                                               #
                data=float(data_pre)                                            #
            else:                                                               #
                data=data_pre                                                   #
            return_list.append(data)                                            #
    if include_randomized_match_data:                                           #
        return return_list, totals, matches                                     #
    else:                                                                       #
        return return_list                                                      #
def sum_active_chanels(string):                                                 #
    return_val=0                                                                #
    for val in string:                                                          #
        assert val=='1' or val=='0'                                             #
        if val == '1':                                                          #
            return_val+=1                                                       #
    return return_val                                                           #
def simple_std(data_like, ddof):                                                #
    processed_data=[]                                                           #
    for row in data_like:                                                       #
        try:                                                                    #
            processed_data.append(float(row))                                   #
        except:                                                                 #
            pass                                                                #
    if len(processed_data) > 1:                                                 #
        return str(stdev(processed_data,ddof=ddof))                             #
    else:                                                                       #
        return 'N/A'                                                            #
def simple_count(data_like):                                                    #
    data_count=0                                                                #
    for row in data_like:                                                       #
        if row != 'N/A':                                                        #
            data_count+=1                                                       #
    return data_count                                                           #
def make_all_colors_order():                                                    #
    tuples=list(itertools.product(['1','0'], repeat=5))[::-1]                   #
    output_list=[]                                                              #
    for row in tuples:                                                          #
        output_list.append(''.join(row))                                        #
    return output_list                                                          #
def count_data(list_of_colors):                                                 #
    number_of_nonblanks=0                                                       #
    for color in list_of_colors:                                                #
        if color != '':                                                         #
            number_of_nonblanks+=1                                              #
    return number_of_nonblanks                                                  #
def count_colors(list_of_colors):                                               #
    number_of_colors=0                                                          #
    for color in list_of_colors:                                                #
        if set(color) != {'0'}:                                                 #
            number_of_colors+=1                                                 #
    return number_of_colors                                                     #
def both_have_data(neuron_1_colors,neuron_2_colors,rep_neurons=False):          #
    if not rep_neurons:                                                         #
        if count_data(neuron_1_colors) > 0 and count_data(neuron_2_colors) > 0: #
            return True                                                         #
        else:                                                                   #
            return False                                                        #
    else:                                                                       #
        if count_data(neuron_1_colors) > 1:                                     #
            return True                                                         #
        else:                                                                   #
            return False                                                        #
def either_have_color(neuron_1_colors,neuron_2_colors,rep_neurons=False):       #
    if not rep_neurons:                                                         #
        if (count_colors(neuron_1_colors) > 0 or                                #
            count_colors(neuron_2_colors) > 0):                                 #
            return True                                                         #
        else:                                                                   #
            return False                                                        #
    else:                                                                       #
        if count_colors(neuron_1_colors) > 0:                                   #
            return True                                                         #
        else:                                                                   #
            return False                                                        #
def both_have_color(neuron_1_colors,neuron_2_colors,rep_neurons=False):         #
    if not rep_neurons:                                                         #
        if (count_colors(neuron_1_colors) > 0 and                               #
            count_colors(neuron_2_colors) > 0):                                 #
            return True                                                         #
        else:                                                                   #
            return False                                                        #
    else:                                                                       #
        if count_colors(neuron_1_colors) > 1:                                   #
            return True                                                         #
        else:                                                                   #
            return False                                                        #
def both_have_match(neuron_1_colors,neuron_2_colors,rep_neurons=False):         #
    if (len(find_matches(neuron_1_colors,                                       #
                        neuron_2_colors,                                        #
                        rep_neurons=rep_neurons)) != 0):                        #
        return True                                                             #
    else:                                                                       #
        return False                                                            #
def find_matches(neuron_1_colors,neuron_2_colors,rep_neurons=False):            #
    if not rep_neurons:                                                         #
        matches=list(set(neuron_1_colors).intersection(neuron_2_colors))        #
    else:                                                                       #
        matches=list(                                                           #
                set(                                                            #
                [x for x in neuron_1_colors if neuron_1_colors.count(x) > 1]))  #
    if '' in matches:                                                           #
        matches.remove('')                                                      #
    if '00000' in matches:                                                      #
        matches.remove('00000')                                                 #
    return matches                                                              #
def randomize_objects_in_table(input_table):                                    #
    ordered_input_list=[]                                                       #
    for row in input_table:                                                     #
        for col in row:                                                         #
            if random_only_colors:                                              #
                if col != '' and col != '00000':                                #
                    ordered_input_list.append(col)                              #
            else:                                                               #
                ordered_input_list.append(col)                                  #
    shuffle(ordered_input_list)                                                 #
    output_table=[]                                                             #
    itter=0                                                                     #
    for row in input_table:                                                     #
        next_line=[]                                                            #
        for col in row:                                                         #
            if random_only_colors:                                              #
                if col != '' and col != '00000':                                #
                    next_line.append(ordered_input_list[itter])                 #
                    itter+=1                                                    #
                else:                                                           #
                    next_line.append(col)                                       #
            else:                                                               #
                next_line.append(ordered_input_list[itter])                     #
                itter+=1                                                        #
        output_table.append(next_line)                                          #
    return output_table                                                         #
def simp_color_match_prob(colors,color_match_prob_dic):                         #
    return_val=0                                                                #
    for color in colors:                                                        #
        if color_match_prob_dic[color] > return_val:                            #
            return_val=color_match_prob_dic[color]                              #
    return return_val                                                           #
def matching_colors_no_rep(line,coords_1,coords_2):                             #
    sublist_1=[]                                                                #
    sublist_2=[]                                                                #
    for coord in coords_1:                                                      #
        if line[coord] != '' and line[coord] != '00000':                        #
            sublist_1.append(line[coord])                                       #
    for coord in coords_2:                                                      #
        if line[coord] != '' and line[coord] != '00000':                        #
            sublist_2.append(line[coord])                                       #
    if set(sublist_1).isdisjoint(sublist_2):                                    #
        return 0                                                                #
    else:                                                                       #
        return 1                                                                #
def not_matching_colors_no_rep(line,coords_1,coords_2):                         #
    sublist_1=[]                                                                #
    sublist_2=[]                                                                #
    sublist_1_pos=False                                                         #
    sublist_2_pos=False                                                         #
    for coord in coords_1:                                                      #
        if line[coord] != '':                                                   #
            sublist_1_pos=True                                                  #
            if line[coord] != '00000':                                          #
                sublist_1.append(line[coord])                                   #
    for coord in coords_2:                                                      #
        if line[coord] != '':                                                   #
            sublist_2_pos=True                                                  #
            if line[coord] != '00000':                                          #
                sublist_2.append(line[coord])                                   #
    if ((sublist_1 != [] or sublist_2 != []) and                                #
        (sublist_1_pos and sublist_2_pos)):                                     #
        return 1                                                                #
    else:                                                                       #
        return 0                                                                #
def matching_colors_prob_no_rep(line,coords_1,coords_2,color_match_prob_dic):   #
    sublist_1=[]                                                                #
    sublist_2=[]                                                                #
    for coord in coords_1:                                                      #
        if line[coord] != '' and line[coord] != '00000':                        #
            sublist_1.append(line[coord])                                       #
    for coord in coords_2:                                                      #
        if line[coord] != '' and line[coord] != '00000':                        #
            sublist_2.append(line[coord])                                       #
    matches=list(set(sublist_1) & set(sublist_2))                               #
    if len(matches) == 0:                                                       #
        return 1                                                                #
    else:                                                                       #
        return simp_color_match_prob(matches,color_match_prob_dic)              #
def matching_colors_rep(line,coords):                                           #
    sublist=[]                                                                  #
    for coord in coords:                                                        #
        if line[coord]!= '' and line[coord] != '00000':                         #
            sublist.append(line[coord])                                         #
    if any(sublist.count(x) > 1 for x in sublist):                              #
        return 1                                                                #
    else:                                                                       #
        return 0                                                                #
def not_matching_colors_rep(line,coords):                                       #
    sublist=[]                                                                  #
    sublist_with_blanks=[]                                                      #
    for coord in coords:                                                        #
        if line[coord]!= '':                                                    #
            sublist_with_blanks.append(line[coord])                             #
            if line[coord] != '00000':                                          #
                sublist.append(line[coord])                                     #
    if len(sublist)>0 and len(sublist_with_blanks)>1:                           #
        return 1                                                                #
    else:                                                                       #
        return 0                                                                #
def matching_colors_prob_rep(line,coords,color_match_prob_dic):                 #
    sublist=[]                                                                  #
    for coord in coords:                                                        #
        if line[coord]!= '' and line[coord] != '00000':                         #
            sublist.append(line[coord])                                         #
    matches=[item for item,                                                     #
             count in collections.Counter(sublist).items() if count > 1]        #
    if len(matches)==0:                                                         #
        return 1                                                                #
    else:                                                                       #
        return simp_color_match_prob(matches,color_match_prob_dic)              #
def get_test_stat(input_table,coords_1,coords_2,color_match_prob_dic):          #
    if coords_1 == coords_2:                                                    #
        same_neuron = True                                                      #
    else:                                                                       #
        same_neuron = False                                                     #
    if include_randomized_match_data:                                           #
        matches=0                                                               #
        both_neurons_observed_one_with_color=0                                  #
    if test_stat == 'number_of_matches':                                        #
        matches=0                                                               #
        for row in input_table:                                                 #
            if same_neuron:                                                     #
                matches+=matching_colors_rep(row,coords_1)                      #
                if include_randomized_match_data:                               #
                    both_neurons_observed_one_with_color+=(                     #
                                          not_matching_colors_rep(row,coords_1))#
            else:                                                               #
                matches+=matching_colors_no_rep(row,coords_1,coords_2)          #
                if include_randomized_match_data:                               #
                    both_neurons_observed_one_with_color+=(                     #
                              not_matching_colors_no_rep(row,coords_1,coords_2))#
        test_stat_value = matches                                               #
    elif test_stat == 'percent_match':                                          #
        test_stat_value='N/A'                                                   #
        matches=0                                                               #
        both_neurons_observed_one_with_color=0                                  #
        for row in input_table:                                                 #
            if same_neuron:                                                     #
                matches+=matching_colors_rep(row,coords_1)                      #
                both_neurons_observed_one_with_color+=(                         #
                                          not_matching_colors_rep(row,coords_1))#
            else:                                                               #
                matches+=matching_colors_no_rep(row,coords_1,coords_2)          #
                both_neurons_observed_one_with_color+=(                         #
                              not_matching_colors_no_rep(row,coords_1,coords_2))#
        if float(both_neurons_observed_one_with_color) != float(0):             #
            test_stat_value=str(float(matches)/                                 #
                              float(both_neurons_observed_one_with_color))      #
        else:                                                                   #
            test_stat_value='N/A'                                               #
    elif test_stat == ('sort_of_probibility_of_matches' or                      #
                       'sort_of_prob_of_matches_times_misses'):                 #
        prob=1                                                                  #
        mult=1                                                                  #
        for row in input_table:                                                 #
            if same_neuron:                                                     #
                prob=prob*matching_colors_prob_rep(row,coords_1,                #
                                                   color_match_prob_dic)        #
                if test_stat == 'sort_of_prob_of_matches_times_misses':         #
                    mult+=(not_matching_colors_rep(row,coords_1)-               #
                           matching_colors_rep(row,coords_1))                   #
                if include_randomized_match_data:                               #
                    matches+=matching_colors_rep(row,coords_1)                  #
                    both_neurons_observed_one_with_color+=(                     #
                                          not_matching_colors_rep(row,coords_1))#
            else:                                                               #
                prob=prob*matching_colors_prob_no_rep(row,                      #
                                                      coords_1,                 #
                                                      coords_2,                 #
                                                      color_match_prob_dic)     #
                if test_stat == 'sort_of_prob_of_matches_times_misses':         #
                    mult+=(not_matching_colors_no_rep(row,coords_1,coords_2)-   #
                           matching_colors_no_rep(row,coords_1,coords_2))       #
                if include_randomized_match_data:                               #
                    matches+=matching_colors_no_rep(row,coords_1,coords_2)      #
                    both_neurons_observed_one_with_color+=(                     #
                              not_matching_colors_no_rep(row,coords_1,coords_2))#
        assert mult > 0                                                         #
        test_stat_value=prob*mult                                               #
    if include_randomized_match_data:                                           #
        return test_stat_value, both_neurons_observed_one_with_color, matches   #
    else:                                                                       #
        return test_stat_value                                                  #
#################################################################################
#                    Defineing my data collection functions                     #
#################################################################################
import csv, math, os                                                            #
cwd=os.getcwd()                                                                 #
if '\\' in cwd:                                                                 #
    slash='\\'                                                                  #
else:                                                                           #
    slash='/'                                                                   #
def load_pval_dic(path):                                                        #
    json_string=open(path).read()                                               #
    simp_dic = json.loads(json_string)                                          #
    pval_dic={}                                                                 #
    pval_dic['itterations']=simp_dic['itterations']                             #
    for key in simp_dic:                                                        #
        if key != 'itterations':                                                #
            pval_dic[str(key)]={}                                               #
            for data_key in simp_dic[key]:                                      #
                if (str(data_key) == 'max_value' or                             #
                    str(data_key) == 'min_value' or                             #
                    str(data_key) == 'totals' or                                #
                    str(data_key) == 'matches'):                                #
                    pval_dic[str(key)][str(data_key)]=simp_dic[key][data_key]   #
                    if str(data_key) == 'totals' or str(data_key) == 'matches': #
                        pval_dic[str(key)][str(data_key)]=(                     #
                                                   int(simp_dic[key][data_key]))#
                else:                                                           #
                    if data_key == 'N/A':                                       #
                        pval_dic[str(key)][str(data_key)]=(                     #
                                                        simp_dic[key][data_key])#
                    else:                                                       #
                        pval_dic[str(key)][float(data_key)]=(                   #
                                                        simp_dic[key][data_key])#
    return pval_dic                                                             #
def find_col_ID(table,search_header):                                           #
    header=table[0]                                                             #
    for col_ID in range(0,len(header)):                                         #
        if header[col_ID] ==  search_header:                                    #
            return col_ID                                                       #
    print('We could not find your column')                                      #
    sys.exit(1)                                                                 #
def the_opener(file_path):                                                      #
    if '.csv' in file_path:                                                     #
        return_list=[]                                                          #
        with open(file_path, 'rb') as csvfile:                                  #
            for row in csv.reader(csvfile, dialect='excel'):                    #
                return_list.append(row)                                         #
    elif '.tsv' in file_path or '.txt' in file_path:                            #
        return_list=[]                                                          #
        fixed_file='\n'.join(open(file_path).read().split('\r'))                #
        for row in fixed_file.split('\n'):                                      #
            if row != '':                                                       #
                return_list.append(row.split('\t'))                             #
    return return_list                                                          #
def the_saver(output_path, list_of_lists):                                      #
    with open(output_path,'wb') as csvfile:                                     #
        for row in list_of_lists:                                               #
            csv.writer(csvfile, dialect='excel').writerow(row)                  #
def is_number(value):                                                           #
    try:                                                                        #
        if not math.isnan(float(value)):                                        #
            return True                                                         #
        else:                                                                   #
            return False                                                        #
    except:                                                                     #
        return False                                                            #
def only_numb_keys(dic):                                                        #
    return_list=[]                                                              #
    for key in dic:                                                             #
        if is_number(key):                                                      #
            return_list.append(key)                                             #
    return_list.sort()                                                          #
    return return_list                                                          #
def exists(directory):                                                          #
    try:                                                                        #
        os.stat(directory)                                                      #
        return True                                                             #
    except:                                                                     #
        return False                                                            #
def mkdir(directory):                                                           # This function makes a directory if it does not exist
    if not exists(directory):                                                   #
        if directory[-1]=='\\' or directory[-1]=='/':                           #
            directory_edit=directory[0:-1]                                      #
        else:                                                                   #
            directory_edit=directory                                            #
        make_directory_list=[]                                                  #
        while not exists(directory_edit):                                       #
            make_directory_list.append(directory_edit)                          #
            directory_edit=slash.join(directory_edit.split(slash)[0:-1])        #
        make_directory_list=make_directory_list[::-1]                           #
        for row in make_directory_list:                                         #
            os.mkdir(row)                                                       #
#################################################################################
#                    generating p-value dictionaries functions                  #
#################################################################################
Triangle_table_denomanator=("Total number of times both neurons were observed "+#
                            "and at least one had a color")                     #
def make_test_stat_dic(data_obj,                                                #
                       number_of_p_val_itterations,                             #
                       neuron_layer,                                            #
                       neurons_dic,                                             #
                       neuron_pairs_name_dic,                                   #
                       return_dic):                                             #
    color_match_prob_dic=data_obj.color_match_prob_dic                          #
    for itter in range(0,number_of_p_val_itterations):                          #
        working_data=data_obj.make_random_table()                               #
        for neuron_set_name in neuron_pairs_name_dic[neuron_layer]:             #
            neuron_1=neuron_set_name.split(' to ')[0]                           #
            neuron_2=neuron_set_name.split(' to ')[1]                           #
            coords_1=neurons_dic[neuron_layer][neuron_1].col_IDs                #
            coords_2=neurons_dic[neuron_layer][neuron_2].col_IDs                #
            if include_randomized_match_data:                                   #
                test_stat,totals,matches=get_test_stat(working_data,            #
                                                       coords_1,                #
                                                       coords_2,                #
                                                       color_match_prob_dic)    #
            else:                                                               #
                test_stat=get_test_stat(working_data,                           #
                                        coords_1,                               #
                                        coords_2,                               #
                                        color_match_prob_dic)                   #
            if neuron_set_name not in return_dic:                               #
                return_dic[neuron_set_name]=''                                  #
            if include_randomized_match_data:                                   #
                return_dic[neuron_set_name]+=(str(test_stat)+'|||'+             #
                                              str(totals)+'|||'+                #
                                              str(matches)+',')                 #
            else:                                                               #
                return_dic[neuron_set_name]+=str(test_stat)+','                 #
#################################################################################
#                       simple to use processing functions                      #
#################################################################################
def make_simple_outputs(data_object,path):                                      #
    if path[-1] != slash:                                                       #
        path+=slash                                                             #
    color_number_stat_output_path=path+'color_number_stats'                     #
    color_stat_path=path+'color_stats'                                          #
    group_color_stat_path=path+'group_color_stats'                              #
    the_saver(color_number_stat_output_path+'.csv',                             #
              data_object.make_color_number_stats())                            #
    the_saver(color_stat_path+'.csv',                                           #
              data_object.make_color_stats())                                   #
    the_saver(group_color_stat_path+'.csv',                                     #
              data_object.make_group_color_stats())                             #
    if not isinstance(data_object,larva):                                       #
        if include_by_the_larva_outputs:                                        #
            the_saver(color_number_stat_output_path+'_by_the_larva.csv',        #
                    data_object.make_color_number_stats_larva())                #
            the_saver(color_stat_path+'_by_the_larva.csv',                      #
                    data_object.make_color_stats_larva())                       #
            the_saver(group_color_stat_path+'_by_the_larva.csv',                #
                    data_object.make_group_color_stats_larva())                 #
def make_all_outputs(data_object,                                               #
                      path,                                                     #
                      flat_outputs_dic,                                         #
                      include_in_flaten=False,                                  #
                      include_larva_data=False,                                 #
                      include_pvalues=False):                                   #
    if path[-1] != slash:                                                       #
        path+=slash                                                             #
    make_simple_outputs(data_object,path)                                       #
    for neuron_layer in neurons_dic.keys():                                     #
        if neuron_layer not in exclude_color_layers:                            #
            if neuron_layer not in flat_outputs_dic and include_in_flaten:      #
                flat_outputs_dic[neuron_layer]={}                               #
            current_path=path+neuron_layer+slash                                #
            mkdir(current_path)                                                 #
            error_bar_dic='N/A'                                                 #
            if include_larva_data:                                              #
                error_bar_dic=make_by_the_colors_output_larva(data_object,      #
                                                              neuron_layer,     #
                                                              current_path)     #
            make_by_the_colors_output(data_object,                              #
                                    neuron_layer,                               #
                                    current_path,                               #
                                    flat_outputs_dic,                           #
                                    include_in_flaten=include_in_flaten,        #
                                    include_pvalues=include_pvalues,            #
                                    error_bar_dic=error_bar_dic)                #
def make_by_the_colors_output(data_object,                                      #
                              neuron_layer,                                     #
                              path,                                             #
                              flat_outputs_dic,                                 #
                              include_in_flaten=False,                          #
                              include_pvalues=False,                            #
                              error_bar_dic='N/A'):                             #
    assert isinstance(data_object, basefunctions)                               #
    working_layer_dic=neurons_dic[neuron_layer]                                 #
    list_of_sets_already_done=[]                                                #
    output_table=''                                                             #
    list_of_neurons=key_file_order_dic[neuron_layer]                            #
    for neuron_1_ID in list_of_neurons:                                         #
        neuron_1=working_layer_dic[neuron_1_ID]                                 #
        for neuron_2_ID in list_of_neurons:                                     #
            neuron_2=working_layer_dic[neuron_2_ID]                             #
            working_set=set([neuron_1_ID,neuron_2_ID])                          #
            if working_set not in list_of_sets_already_done:                    #
                list_of_sets_already_done.append(working_set)                   #
                if output_table == '':                                          #
                    output_table=data_object.make_match_report(neuron_1,        #
                                                              neuron_2,         #
                                                            include_header=True,#
                                                include_pvalues=include_pvalues,#
                                                             pval_dic_path=path)#
                else:                                                           #
                    output_table.append(data_object.make_match_report(neuron_1, #
                                                                      neuron_2, #
                                               include_pvalues=include_pvalues, #
                                                            pval_dic_path=path))#
    if include_pvalues:                                                         #
        p_val_column=output_table[0].index('Relatedness test P-Value')          #
        p_val_list=[]                                                           #
        for row in output_table[1:]:                                            #
            p_val_list.append(float(row[p_val_column]))                         #
        corrected_ps=multipletests(p_val_list,                                  #
                                   alpha=alpha_level,                           #
                                   method=Multiple_hypothesis_method)[1]        #
        output_table[0]=(output_table[0][0:p_val_column+1]+                     #
                         ['FDR']+                                               #
                         output_table[0][p_val_column+1:])                      #
        for row_ID in range(1,len(output_table)):                               #
            output_table[row_ID]=(output_table[row_ID][0:p_val_column+1]+       #
                                  [corrected_ps[row_ID-1]]+                     #
                                  output_table[row_ID][p_val_column+1:])        #
    the_saver(path+'by_the_colors.csv',output_table)                            #
    p_val_col='N/A'                                                             #
    pair_name_col='N/A'                                                         #
    both_obs_col='N/A'                                                          #
    either_color_col='N/A'                                                      #
    both_color_col='N/A'                                                        #
    match_col='N/A'                                                             #
    for col_head_ID in range(0,len(output_table[0])):                           #
        if (include_pvalues and                                                 #
            output_table[0][col_head_ID] == 'Relatedness test P-Value'):        #
            p_val_col=col_head_ID                                               #
        if output_table[0][col_head_ID] == 'Simple Match Name':                 #
            pair_name_col=col_head_ID                                           #
        if output_table[0][col_head_ID] ==('Total number of times both neurons'+#
                                           ' were observed'):                   #
            both_obs_col=col_head_ID                                            #
            if Triangle_table_denomanator == output_table[0][col_head_ID]:      #
                dem_col=col_head_ID                                             #
        if output_table[0][col_head_ID] ==('Total number of times both neurons'+#
                                 ' were observed and at least one had a color'):#
            either_color_col=col_head_ID                                        #
            if Triangle_table_denomanator == output_table[0][col_head_ID]:      #
                dem_col=col_head_ID                                             #
        if output_table[0][col_head_ID] ==('Total number of times both neurons'+#
                                           ' had color'):                       #
            both_color_col=col_head_ID                                          #
            if Triangle_table_denomanator == output_table[0][col_head_ID]:      #
                dem_col=col_head_ID                                             #
        if output_table[0][col_head_ID] ==('Total number of times both neurons'+#
                                           ' had the same color (matches)'):    #
            match_col=col_head_ID                                               #
            if Triangle_table_denomanator == output_table[0][col_head_ID]:      #
                dem_col=col_head_ID                                             #
        if output_table[0][col_head_ID] == 'FDR':                               #
            p_val_col_corrected=col_head_ID                                     #
    p_val_pair_dic={}                                                           #
    both_obs_pair_dic={}                                                        #
    either_color_pair_dic={}                                                    #
    both_color_pair_dic={}                                                      #
    match_col_pair_dic={}                                                       #
    dem_col_pair_dic={}                                                         #
    corrected_p_val_dic={}                                                      #
    for row in output_table[1:]:                                                #
        pair_name=row[pair_name_col]                                            #
        both_obs=row[both_obs_col]                                              #
        both_obs_pair_dic[pair_name]=both_obs                                   #
        either_color=row[either_color_col]                                      #
        either_color_pair_dic[pair_name]=either_color                           #
        both_color=row[both_color_col]                                          #
        both_color_pair_dic[pair_name]=both_color                               #
        matches=row[match_col]                                                  #
        match_col_pair_dic[pair_name]=matches                                   #
        dem_value=row[dem_col]                                                  #
        dem_col_pair_dic[pair_name]=dem_value                                   #
        if include_pvalues:                                                     #
            p_val=row[p_val_col]                                                #
            p_val_pair_dic[pair_name]=p_val                                     #
            corrected_p_val_dic[pair_name]=row[p_val_col_corrected]             #
        neuron_A_name=pair_name.split(' to ')[0]                                #
    neuron_num=len(list_of_neurons)                                             #
    assert neuron_num*(neuron_num-1)/2+neuron_num+1 == len(output_table)        # makes sure all the neurons are accounted for
    if include_pvalues:                                                         #
        p_val_outputs=[['']+list_of_neurons]                                    #
        p_val_corected_outputs=[['']+list_of_neurons]                           #
    both_obs_outputs=[['']+list_of_neurons]                                     #
    either_color_outputs=[['']+list_of_neurons]                                 #
    both_color_outputs=[['']+list_of_neurons]                                   #
    match_outputs=[['']+list_of_neurons]                                        #
    triangle_table_frac=[['']+list_of_neurons]                                  #
    trinagle_table_percent=[['']+list_of_neurons]                               #
    for neuron_A_name in list_of_neurons:                                       #
        if include_pvalues:                                                     #
            next_line_p=[neuron_A_name]                                         #
        next_line_both_obs=[neuron_A_name]                                      #
        next_line_either_color=[neuron_A_name]                                  #
        next_line_both_color=[neuron_A_name]                                    #
        next_line_match=[neuron_A_name]                                         #
        next_line_triangle_table_frac=[neuron_A_name]                           #
        next_line_triangle_table_percent=[neuron_A_name]                        #
        next_line_corrected_p=[neuron_A_name]                                   #
        for neuron_B_name in list_of_neurons:                                   #
            neuron_match_key=neuron_A_name+' to '+neuron_B_name                 #
            if neuron_match_key in both_obs_pair_dic:                           #
                if include_pvalues:                                             #
                    next_line_p.append(p_val_pair_dic[neuron_match_key])        #
                    next_line_corrected_p.append(corrected_p_val_dic[           #
                                                              neuron_match_key])#
                next_line_both_obs.append(both_obs_pair_dic[neuron_match_key])  #
                next_line_either_color.append(                                  #
                                              either_color_pair_dic[            #
                                              neuron_match_key])                #
                next_line_both_color.append(                                    #
                                              both_color_pair_dic[              #
                                              neuron_match_key])                #
                next_line_match.append(                                         #
                                              match_col_pair_dic[               #
                                              neuron_match_key])                #
                next_line_triangle_table_frac.append('"'+str(                   #
                                                    match_col_pair_dic[         #
                                                    neuron_match_key])+"/"+     #
                                                    str(                        #
                                                    dem_col_pair_dic[           #
                                                    neuron_match_key])+'"')     #
                if str(dem_col_pair_dic[neuron_match_key]) != '0':              #
                    next_line_triangle_table_percent.append(str(                #
                                                            float(              #
                                                            match_col_pair_dic[ #
                                                            neuron_match_key])/ #
                                                            float(              #
                                                            dem_col_pair_dic[   #
                                                            neuron_match_key])))#
                else:                                                           #
                    next_line_triangle_table_percent.append("N/A")              #
            else:                                                               #
                if include_pvalues:                                             #
                    next_line_p.append('')                                      #
                    next_line_corrected_p.append('')                            #
                next_line_both_obs.append('')                                   #
                next_line_either_color.append('')                               #
                next_line_both_color.append('')                                 #
                next_line_match.append('')                                      #
                next_line_triangle_table_frac.append('')                        #
                next_line_triangle_table_percent.append('')                     #
        if include_pvalues:                                                     #
            p_val_outputs.append(next_line_p)                                   #
            p_val_corected_outputs.append(next_line_corrected_p)                #
        both_obs_outputs.append(next_line_both_obs)                             #
        either_color_outputs.append(next_line_either_color)                     #
        both_color_outputs.append(next_line_both_color)                         #
        match_outputs.append(next_line_match)                                   #
        triangle_table_frac.append(next_line_triangle_table_frac)               #
        trinagle_table_percent.append(next_line_triangle_table_percent)         #
    if include_pvalues:                                                         #
        the_saver(path+'pval_matrix.csv',p_val_outputs)                         #
        the_saver(path+'FDR_matrix.csv',p_val_corected_outputs)                 #
    the_saver(path+'Match_fraction_matrix.csv',triangle_table_frac)             #
    the_saver(path+'Match_precent_matrix.csv',trinagle_table_percent)           #
    data_name=data_object.name                                                  #
    if include_in_flaten:                                                       #
        if include_pvalues:                                                     #
            if 'p_val_matrix' not in flat_outputs_dic[neuron_layer]:            #
                flat_outputs_dic[neuron_layer]['p_val_matrix']={}               #
            flat_outputs_dic[neuron_layer]['p_val_matrix'][data_name]=(         #
                                                                  p_val_outputs)#
        if 'both_obs_matrix' not in flat_outputs_dic[neuron_layer]:             #
            flat_outputs_dic[neuron_layer]['both_obs_matrix']={}                #
        flat_outputs_dic[neuron_layer]['both_obs_matrix'][data_name]=(          #
                                                               both_obs_outputs)#
        if 'either_color_matrix' not in flat_outputs_dic[neuron_layer]:         #
            flat_outputs_dic[neuron_layer]['either_color_matrix']={}            #
        flat_outputs_dic[neuron_layer]['either_color_matrix'][data_name]=(      #
                                                           either_color_outputs)#
        if 'both_color_matrix' not in flat_outputs_dic[neuron_layer]:           #
            flat_outputs_dic[neuron_layer]['both_color_matrix']={}              #
        flat_outputs_dic[neuron_layer]['both_color_matrix'][data_name]=(        #
                                                             both_color_outputs)#
        if 'match_matrix' not in flat_outputs_dic[neuron_layer]:                #
            flat_outputs_dic[neuron_layer]['match_matrix']={}                   #
        flat_outputs_dic[neuron_layer]['match_matrix'][data_name]=match_outputs #
        if 'triangle_table_frac' not in flat_outputs_dic[neuron_layer]:         #
            flat_outputs_dic[neuron_layer]['triangle_table_frac']={}            #
        flat_outputs_dic[neuron_layer]['triangle_table_frac'][data_name]=(      #
                                                            triangle_table_frac)#
        if 'triangle_table_percent' not in flat_outputs_dic[neuron_layer]:      #
            flat_outputs_dic[neuron_layer]['triangle_table_percent']={}         #
        flat_outputs_dic[neuron_layer]['triangle_table_frac'][data_name]=(      #
                                                         trinagle_table_percent)#
def make_by_the_colors_output_larva(data_object,neuron_layer,path):             #
    error_bar_dic={}                                                            #
    if path[-1] != slash:                                                       #
        path+=slash                                                             #
    working_path=path+'By the larva'+slash                                      #
    mkdir(working_path)                                                         #
    assert isinstance(data_object, basefunctions)                               #
    working_layer_dic=neurons_dic[neuron_layer]                                 #
    list_of_sets_already_done=[]                                                #
    list_of_neurons=key_file_order_dic[neuron_layer]                            #
    for neuron_1_ID in list_of_neurons:                                         #
        neuron_1=working_layer_dic[neuron_1_ID]                                 #
        for neuron_2_ID in list_of_neurons:                                     #
            neuron_2=working_layer_dic[neuron_2_ID]                             #
            working_set=set([neuron_1_ID,neuron_2_ID])                          #
            if working_set not in list_of_sets_already_done:                    #
                list_of_sets_already_done.append(working_set)                   #
                simple_name_pre=neuron_1.identifier + ' to '+neuron_2.identifier#
                simple_name='-'.join(simple_name_pre.split('/'))                #
                working_output_path=working_path+simple_name+'.csv'             #
                output=data_object.make_match_report_larva(neuron_1,neuron_2)   #
                the_saver(working_output_path,output)                           #
                error_bar_column=find_col_ID(output,'Fractional match')         #
                error_bar_dic[simple_name]=[]                                   #
                for row in output[1:]:                                          #
                    error_bar_dic[simple_name].append(row[error_bar_column])    #
    return error_bar_dic                                                        #
def check_if_in(simple_hemisegment_str,included_segments):                      #
    for included_segment in included_segments:                                  #
        if included_segment in simple_hemisegment_str:                          #
            return True                                                         #
    return False                                                                #
#################################################################################
#                        multithreaded p_val functions                          #
################################################################################# 
import multiprocessing, os, sys                                                 #
from time import sleep                                                          #
def _getThreads():                                                              #
    " Returns the number of available threads on a posix/win based system "     #
    if sys.platform == 'win32':                                                 #
        return (int)(os.environ['NUMBER_OF_PROCESSORS'])                        #
    else:                                                                       #
        return (int)(os.popen('grep -c cores /proc/cpuinfo').read())            #
max_threads=int(round(_getThreads()*Thread_multiplyer,0))+Thread_offset         #
worker_dic={}                                                                   #
finished_workers_dic={}                                                         #
for thread in range(0,max_threads):                                             #
    worker_dic[thread]=''                                                       #
    finished_workers_dic[thread]=[True,'']                                      #
def jobs_done(jobs_dic):                                                        #
    for job_ID in jobs_dic:                                                     #
        if jobs_dic[job_ID] != 'Done':                                          #
            return False                                                        #
    return True                                                                 #
def add_two_simple_dics(simp_dic_1,simp_dic_2):                                 #
    if simp_dic_1 == {}:                                                        #
        assert simp_dic_2 != {}                                                 #
        return simp_dic_2                                                       #
    elif simp_dic_2 == {}:                                                      #
        return simp_dic_1                                                       #
    itter_1=simp_dic_1['itterations']                                           #
    itter_2=simp_dic_2['itterations']                                           #
    assert len(simp_dic_1)==len(simp_dic_2)                                     #
    new_itter=itter_1+itter_2                                                   #
    return_dic={}                                                               #
    return_dic['itterations']=new_itter                                         #
    for neuron_match_name in simp_dic_1:                                        #
        if neuron_match_name != 'itterations':                                  #
            return_dic[neuron_match_name]={}                                    #
            data_from_1=simp_dic_1[neuron_match_name]                           #
            data_from_2=simp_dic_2[neuron_match_name]                           #
            posible_values=list(set(data_from_1.keys()+data_from_2.keys()))     #
            posible_values.remove('max_value')                                  #
            posible_values.remove('min_value')                                  #
            if include_randomized_match_data:                                   #
                posible_values.remove('totals')                                 #
                posible_values.remove('matches')                                #
            for posible_value in posible_values:                                #
                count_of_value=0                                                #
                if posible_value in simp_dic_1[neuron_match_name]:              #
                    count_of_value+=simp_dic_1[neuron_match_name][posible_value]#
                if posible_value in simp_dic_2[neuron_match_name]:              #
                    count_of_value+=simp_dic_2[neuron_match_name][posible_value]#
                return_dic[neuron_match_name][posible_value]=count_of_value     #
            return_dic[neuron_match_name]['max_value']=max(posible_values)      #
            return_dic[neuron_match_name]['min_value']=min(posible_values)      #
            if include_randomized_match_data:                                   #
                return_dic[neuron_match_name]['totals']=(                       #
                              simp_dic_1[neuron_match_name]['totals']+          #
                              simp_dic_2[neuron_match_name]['totals'])          #
                return_dic[neuron_match_name]['matches']=(                      #
                              simp_dic_1[neuron_match_name]['matches']+         #
                              simp_dic_2[neuron_match_name]['matches'])         #
    return return_dic                                                           #
def finalize_dic(master_dic):                                                   #
    if test_stat=='number_of_matches' or test_stat=='percent_match':            #
        big_or_small_vals_are_good='big'                                        #
    elif (test_stat=='sort_of_probibility_of_matches' or                        #
          test_stat=='sort_of_prob_of_matches_times_misses'):                   #
        big_or_small_vals_are_good='small'                                      #
    else:                                                                       #
        print('I do not understand your test stat '+test_stat)                  #
    final_dic={}                                                                #
    itterations=master_dic['itterations']                                       #
    final_dic['itterations']=itterations                                        #
    for key in master_dic:                                                      #
        if key != 'itterations':                                                #
            final_dic[key]={}                                                   #
            working_data=master_dic[key]                                        #
            posible_values=working_data.keys()                                  #
            posible_values.remove('max_value')                                  #
            posible_values.remove('min_value')                                  #
            if include_randomized_match_data:                                   #
                posible_values.remove('totals')                                 #
                posible_values.remove('matches')                                #
            if big_or_small_vals_are_good == 'small':                           #
                posible_values.sort()                                           #
                final_dic[key]['max_value']=posible_values[-1]                  #
                final_dic[key]['min_value']=posible_values[0]                   #
            else:                                                               #
                posible_values.sort(reverse=True)                               #
                final_dic[key]['max_value']=posible_values[0]                   #
                final_dic[key]['min_value']=posible_values[-1]                  #
            running_total=0                                                     #
            for value in posible_values:                                        #
                if value != 'max_value' and value != 'min_value':               #
                    running_total+=master_dic[key][value]                       #
                    final_dic[key][value]=(float(running_total)/                #
                                           float(itterations))                  #
            assert running_total == itterations                                 #
            if include_randomized_match_data:                                   #
                final_dic[key]['totals']=master_dic[key]['totals']              #
                final_dic[key]['matches']=master_dic[key]['matches']            #
    return final_dic                                                            #
def make_test_stat_dic_multithreaded(data_obj,                                  #
                                     itteration_print,                          #
                                     neuron_layer,                              #
                                     neurons_dic,                               #
                                     neuron_pairs_name_dic,number_of_jobs_pre,  #
                                     pval_dic_path):                            #
    master_dic={}                                                               #
    jobs_dic={}                                                                 #
    return_data_dic={}                                                          #
    temp_output_path=pval_dic_path+'temp_dic.json'                              #
    if os.path.isfile(temp_output_path):                                        #
        master_dic=load_pval_dic(temp_output_path)                              #
        current_itter=master_dic['itterations']                                 #
        assert current_itter % itteration_print == 0                            #
        jobs_done_sucessfully_before=int(current_itter/itteration_print)        #
    else:                                                                       #
        jobs_done_sucessfully_before=0                                          #
    number_of_jobs=number_of_jobs_pre-jobs_done_sucessfully_before              #
    manager = multiprocessing.Manager()                                         #
    for job_ID in range(number_of_jobs):                                        #
        jobs_dic[job_ID]='not started'                                          #
    searching_for_worker=True                                                   #
    queued_job_ID=0                                                             #
    while searching_for_worker:                                                 #
        started_a_job=False                                                     #
        for worker_name in worker_dic:                                          #
            start_a_job_with_this_worker=False                                  #
            if (worker_dic[worker_name] == '' and                               #
                queued_job_ID != number_of_jobs):                               #
                start_a_job_with_this_worker=True                               #
            elif (worker_dic[worker_name] != '' and                             #
                  not worker_dic[worker_name].is_alive()):                      #
                if not finished_workers_dic[worker_name][0]:                    #
                    return_dic_multi = dict(return_data_dic[worker_name])       #
                    simp_dic_multi=simplify_dic(return_dic_multi,               #
                                                neuron_pairs_name_dic,          #
                                                neuron_layer,                   #
                                                itteration_print)               #
                    master_dic=add_two_simple_dics(master_dic,simp_dic_multi)   #
                    json_text_dic=json.dumps(master_dic)                        #
                    open(temp_output_path,'w').write(json_text_dic)             #
                    finished_workers_dic[worker_name][0]=True                   #
                    finished_job_ID=finished_workers_dic[worker_name][1]        #
                    jobs_dic[finished_job_ID]='Done'                            #
                    print('\t\tFinished '+str((finished_job_ID+1+               #
                                               jobs_done_sucessfully_before)*   #
                                               itteration_print)+' iterations.')#
                if not jobs_done(jobs_dic) and queued_job_ID != number_of_jobs: #
                    start_a_job_with_this_worker=True                           #
            if start_a_job_with_this_worker:                                    #
                assert queued_job_ID != number_of_jobs                          #
                return_data_dic[worker_name] = manager.dict()                   #
                return_dic_multi = return_data_dic[worker_name]                 #
                command_tup=(data_obj,                                          #
                             itteration_print,                                  #
                             neuron_layer,                                      #
                             neurons_dic,                                       #
                             neuron_pairs_name_dic,                             #
                             return_dic_multi)                                  #
                worker_dic[worker_name]=multiprocessing.Process(                #
                                                      target=make_test_stat_dic,#
                                                      args=command_tup)         #
                worker_dic[worker_name].start()                                 #
                finished_workers_dic[worker_name][0]=False                      #
                finished_workers_dic[worker_name][1]=queued_job_ID              #
                jobs_dic[queued_job_ID]='Started'                               #
                print('\t\tStarted '+str((queued_job_ID+1+                      #
                                          jobs_done_sucessfully_before)*        #
                                          itteration_print)+' iterations.')     #
                queued_job_ID+=1                                                #
                started_a_job=True                                              #
        if not started_a_job:                                                   #
            if verbose_multithread_worker_search:                               #
                print('\t\twaiting for free workers')                           #
            sleep(multithread_sleep_time)                                       #
            if not jobs_done(jobs_dic):                                         #
                searching_for_worker=True                                       #
            else:                                                               #
                searching_for_worker=False                                      #
    final_dic=finalize_dic(master_dic)                                          #
    os.remove(temp_output_path)                                                 #
    return final_dic                                                            #
def convert_to_P_val(value,list_of_all_values,big_or_small_numbers_are_good):   #
    if big_or_small_numbers_are_good == 'big':                                  #
        count_of_better=sum(i <= value for i in list_of_all_values)             #
    if big_or_small_numbers_are_good == 'small':                                #
        count_of_better=sum(i >= value for i in list_of_all_values)             #
    return float(count_of_better)/float(len(list_of_all_values))                #
def simplify_dic(rand_dic,neuron_pairs_name_dic,neuron_layer,itters):           #
    simp_dic={}                                                                 #
    simp_dic['itterations']=itters                                              #
    for key in rand_dic:                                                        #
        assert key not in simp_dic                                              #
        if include_randomized_match_data:                                       #
            data,totals,matches=float_list_special(rand_dic[key].split(','))    #
        else:                                                                   #
            data=float_list_special(rand_dic[key].split(','))                   #
        simp_dic[key]={}                                                        #
        simp_dic[key]['max_value']=max(data)                                    #
        simp_dic[key]['min_value']=min(data)                                    #
        if include_randomized_match_data:                                       #
            simp_dic[key]['totals']=totals                                      #
            simp_dic[key]['matches']=matches                                    #
        for value in data:                                                      #
            if value not in simp_dic[key]:                                      #
                simp_dic[key][value]=data.count(value)                          #
    return simp_dic                                                             #
#################################################################################
#                      import brainbow data into objects                        #
#################################################################################
if __name__ == "__main__":                                                      #
    included_segment_sets_pre=the_opener(included_hemisegments_path)            #
    included_segment_sets=[]                                                    #
    for row in included_segment_sets_pre:                                       #
        next_line=[]                                                            #
        for col in row:                                                         #
            if col != '':                                                       #
                next_line.append(col)                                           #
        included_segment_sets.append(next_line)                                 #
    for included_segments in included_segment_sets:                             #
        mkdir('.'+slash+'Input')                                                #
        mkdir('.'+slash+'Output '+','.join(included_segments))                  #
        input_file_path='.'+slash+'Input'+slash                                 #
        list_of_input_files=os.listdir(input_file_path)                         #
        if list_of_input_files == []:                                           #
            print('Plese put your input files in the Input folder')             #
            sys.exit(0)                                                         #
        conditions_dic={}                                                       # These dictionarys will serve as look up tables for particular objects
        dates_dic={}                                                            #
        larvae_dic={}                                                           #
        hemisegments_dic={}                                                     #
        for input_file in list_of_input_files:                                  # itterate through all condition input files
            if input_file[-4:]=='.csv':                                         #
                working_input_file_path=input_file_path+input_file              # grabbing the full path of the file I would like to work with
                condition_str=input_file[0:-4]                                  # defining the key that will identify this condition file in the conditions dic
                conditions_dic[condition_str]=condition(condition_str)          # defining my condition object
                working_data=the_opener(working_input_file_path)                # import the data from a CSV to a list of lists
                for row in working_data[1:]:                                    # now I itterate through the current condition file I have open extracting data like date, larva, and segment information
                    simple_date_str=row[0]                                      #
                    simple_larva_str=row[1]                                     #
                    simple_hemisegment_str=row[2]+'_'+row[3]                    #
                    if ((not only_include_some_segments) or                     #
                      (check_if_in(simple_hemisegment_str, included_segments))):#
                        date_str=condition_str+'|||'+simple_date_str            # As I am going to require uneque identifiers for all different types of objects, I concatanate the names as the objects become more and more spicific, starting with concatinating the condition object key with the date object key
                        larva_str=date_str+'|||'+simple_larva_str               #
                        hemisegment_str=larva_str+'|||'+simple_hemisegment_str  #
                        if date_str not in dates_dic:                           # In the next set of lines, I define objects if they have not yet been defined.
                            dates_dic[date_str]=date(date_str,                  #
                                                  conditions_dic[condition_str])#
                        if larva_str not in larvae_dic:                         #
                            larvae_dic[larva_str]=larva(larva_str,              #
                                                        dates_dic[date_str])    #
                        if hemisegment_str not in hemisegments_dic:             #
                            hemisegments_dic[hemisegment_str]=hemisegment(      #
                                                          hemisegment_str,      #
                                                          larvae_dic[larva_str])#
                        else:                                                   #
                            warnings.warn('For some reason you are trying to '  #
                                  'add the same hemisegment object twice, plese'#
                                  ' go back and check your code larva '         #
                                  'object, just so you know.',                  #
                                      SyntaxWarning)                            #
                        neuron_data=row[4:49]                                   #
                        trouble_str=', '.join(hemisegment_str.split('|||')[1:]) #
                        assert len(neuron_data)==45,('check '+input_file+' on ' #
                                  'the line with the hemisegment '+             #
                                  trouble_str+' it '                            #
                                  'looks like there is not data for all 45 '    #
                                  'neurons, check the end of the line. That is '#
                                  'typically where the problems are')           #
                        for col in neuron_data:                                 #
                            if col != '':                                       #
                                assert len(col)==5,('check '+input_file+' on '  #
                             'the line with the hemisegment '+trouble_str+' it '#
                             'looks like one of the cells does not have all '   #
                             '5 color code')                                    #
                                assert (set(col)=={'0', '1'} or                 #
                                    set(col)=={'0'} or                          #
                                    set(col)=={'1'}), (                         #
                                    'check '+input_file+' on the line with the '#
                                    'hemisegment '+trouble_str+' it looks like '#
                                    'the color code contains more than just 1 ' #
                                    'and 0')                                    #
                        hemisegments_dic[hemisegment_str].add_data(neuron_data) #
#################################################################################
#                     inporting key data into neuron objecst                    #
#################################################################################
        key_file=the_opener(path_of_key_file)                                   #
        neurons_dic={}                                                          #
        key_file_order_dic={}                                                   #
        for row in key_file[1:]:                                                #
            for col_ID in range(1,len(key_file[0])):                            #
                neuron_name=row[col_ID]                                         #
                color_layer=key_file[0][col_ID]                                 #
                if color_layer not in key_file_order_dic:                       #
                    key_file_order_dic[color_layer]=[]                          #
                if neuron_name not in key_file_order_dic[color_layer]:          #
                    key_file_order_dic[color_layer].append(neuron_name)         #
                neuron_identifier=int(row[0])-1                                 #
                if color_layer not in neurons_dic:                              #
                    neurons_dic[color_layer]={}                                 #
                if neuron_name not in neurons_dic[color_layer]:                 #
                    neurons_dic[color_layer][neuron_name]=neuron(neuron_name,   #
                                                                color_layer)    #
                neurons_dic[color_layer][neuron_name].add_col_ID(               #
                                                              neuron_identifier)#
        neuron_pairs_name_dic={}                                                #
        for neuron_layer in key_file_order_dic:                                 #
            list_of_neuron_sets=[]                                              #
            neuron_pairs_name_dic[neuron_layer]=[]                              #
            for neuron_1_name in key_file_order_dic[neuron_layer]:              #
                for neuron_2_name in key_file_order_dic[neuron_layer]:          #
                    neuron_set=set([neuron_1_name,neuron_2_name])               #
                    if neuron_set not in list_of_neuron_sets:                   #
                       list_of_neuron_sets.append(neuron_set)                   #
                       simp_set_name = neuron_1_name + ' to '+neuron_2_name     #
                       neuron_pairs_name_dic[neuron_layer].append(simp_set_name)#
#################################################################################
#                                  main code                                    #
#################################################################################
        flat_outputs_dic={}                                                     #
        for condition_identifier in conditions_dic:                             #
            output_path='.'+slash+'Output '+','.join(included_segments)+slash   #
            simple_condition_string=condition_identifier                        #
            print('Processing the '+simple_condition_string+' condition.')      #
            current_path=(output_path+                                          #
                        simple_condition_string+slash+'all_data')               #
            mkdir(current_path)                                                 #
            contisions_object=conditions_dic[condition_identifier]              #
            make_all_outputs(contisions_object,                                 #
                            current_path,                                       #
                            flat_outputs_dic,                                   #
                            include_in_flaten=include_conditions_in_flaten,     #
                            include_larva_data=include_larva_data_conditions,   #
                            include_pvalues=include_pvalue_conditions)          #
            date_objects=contisions_object.dates                                #
            for date_object in date_objects:                                    #
                simple_date_string=date_object.date.split('|||')[-1]            #
                print('\t'+simple_date_string)                                  #
                current_path=(output_path+                                      #
                            simple_condition_string+slash+                      #
                            simple_date_string+slash+'all_data')                #
                mkdir(current_path)                                             #
                make_all_outputs(date_object,                                   #
                                current_path,                                   #
                                flat_outputs_dic,                               #
                                include_in_flaten=include_date_in_flaten,       #
                                include_larva_data=include_larva_data_dates,    #
                                include_pvalues=include_pvalue_dates)           #
                larva_objects=date_object.larvae                                #
                for larva_object in larva_objects:                              #
                    simple_larva_string=larva_object.larva.split('|||')[-1]     #
                    print('\t\t'+simple_larva_string)                           #
                    current_path=(output_path+                                  #
                                simple_condition_string+slash+                  #
                                simple_date_string+slash+                       #
                                simple_larva_string+slash)                      #
                    mkdir(current_path)                                         #
                    make_all_outputs(larva_object,                              #
                                    current_path,                               #
                                    flat_outputs_dic,                           #
                                    include_in_flaten=include_larva_in_flaten,  #
                                    include_larva_data=False,                   #
                                    include_pvalues=include_pvalue_larva)       #
        if len(flat_outputs_dic)!=0:                                            #
             output_path=('.'+slash+                                            #
                         'Output '+                                             #
                        ','.join(included_segments)+                            #
                        slash+                                                  #
                        'Flat_outputs'+                                         #
                        slash)                                                  #
             mkdir(output_path)                                                 #
             for neuron_layer_key in flat_outputs_dic:                          #
                current_path=output_path+neuron_layer_key+slash                 #
                mkdir(current_path)                                             #
                for output_file_name in flat_outputs_dic[neuron_layer_key]:     #
                    output=[]                                                   #
                    current_file_path=current_path+output_file_name+'.csv'      #
                    for data_header_text in flat_outputs_dic[                   #
                                                            neuron_layer_key][  #
                                                            output_file_name]:  #
                        output.append([data_header_text])                       #
                        output+=flat_outputs_dic[                               #
                                                neuron_layer_key][              #
                                                output_file_name][              #
                                                data_header_text]               #
                    the_saver(current_file_path,output)                         #
