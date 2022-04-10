import pandas as pd
import argparse
import data_process


def main(args_main):

    prediction = args_main.prediction
    epi = args_main.epigenomics
    res = args_main.res_dir
    model_dir = args_main.model_dir
    if prediction:
        bedfile = args_main.bed
        df = pd.read_csv(bedfile,sep='\t', header=None)
        df.columns = ['chr', 'pos', 'strand', 'ref', 'alt']
        if epi =='H3K9':
            import Model_histone
            res_histone, ref_marix_histone, alt_martix_histone =data_process.process_data_histone(df,res)


            ref_prediction_histone = Model_histone.main(ref_marix_histone,model_dir)
            alt_prediction_histone = Model_histone.main(alt_martix_histone,model_dir)

            res_histone['H3K9AC_REF_Pred'] = ref_prediction_histone.flatten()
            res_histone['H3K9AC_ALT_Pred'] = alt_prediction_histone.flatten()
            res_histone['DELTA_H3K9AC'] = alt_prediction_histone.flatten() - ref_prediction_histone.flatten()
            res_histone.to_csv('{}/Results_histone.csv'.format(res),sep='\t',index=False)



        if epi =='DNA_Methylation':
            import Model_methy
            res_methy, ref_marix_methy, alt_martix_methy = data_process.process_data_methy(df, res)


            ref_prediction_methy = Model_methy.main(ref_marix_methy,model_dir)
            alt_prediction_methy = Model_methy.main(alt_martix_methy, model_dir)

            res_methy['Methylation_REF_Pred'] = ref_prediction_methy.flatten()
            res_methy['Methylation_ALT_Pred'] = alt_prediction_methy.flatten()
            res_methy['DELTA_Methylation'] = alt_prediction_methy.flatten() - ref_prediction_methy.flatten()
            res_methy.to_csv('{}/Results_methy.csv'.format(res),sep='\t',index=False)


def parse_arguments(parser):

    parser.add_argument('--prediction', dest='prediction', action='store_true', help='Use this option for predict DeepPerVar score')
    parser.add_argument('--epigenomics', type=str, default='H3K9', help='Epigenetics, can be H3K9 or DNA_methylation')
    parser.set_defaults(prediction=False)
    parser.add_argument('--bed', type=str, default='data/snps.bed',
                        help='The Bed file for predicts epigenetics and mutation effects')
    parser.add_argument('--model_dir', type=str, default='/Users/ywang17/workspace/DeepPerVar/models', metavar='<data_directory>',
                        help='The model directory for DeepPerVar')
    parser.add_argument('--res_dir', type=str, default='/Users/ywang17/workspace/DeepPerVar/res', metavar='<data_directory>',
                        help='The data directory for save res')

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='DeepPerVar: a multimodal deep learning framework for functional interpretation of genetic variants in personal genome')
    args_main = parse_arguments(parser)
    main(args_main)