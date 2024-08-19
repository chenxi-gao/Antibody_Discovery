import re
import sys
import Logger
import TimeTracker
import pandas as pd
from Bio import pairwise2
from Bio.Align import substitution_matrices


def main(input_target_species, input_file_path):
    # 验证物种是否为指定列表中的一个
    valid_species = ['fly', 'frog', 'worm', 'mosquito', 'fish']
    verified_target_species = input_target_species.lower()
    if verified_target_species not in valid_species:
        raise ValueError(f"Invalid species '{input_target_species}'. Please choose from {valid_species}.")

    # mode = 'test'
    mode = 'prod'

    time_tracker = TimeTracker.TimeTracker()

    start_time = time_tracker.get_timestamp()

    # database for csv files
    uniprot_human_protein_database_path = 'database/AntibodyFilter_protein_Human.csv'
    uniprot_mouse_protein_database_path = 'database/AntibodyFilter_protein_Mouse.csv'
    uniprot_rat_protein_database_path = 'database/AntibodyFilter_protein_Rat.csv'

    uniprot_protein_database_path = {'human': uniprot_human_protein_database_path,
                                     'mouse': uniprot_mouse_protein_database_path,
                                     'rat': uniprot_rat_protein_database_path}

    human_protein_database_path = 'database/Protein_Human.csv'
    mouse_protein_database_path = 'database/Protein_Mouse.csv'
    rat_protein_database_path = 'database/Protein_Rat.csv'
    fish_protein_database_path = 'database/Protein_Fish.csv'
    fly_protein_database_path = 'database/Protein_Fly.csv'
    frog_protein_database_path = 'database/Protein_Frog.csv'
    mosquito_protein_database_path = 'database/Protein_Mosquito.csv'

    protein_database_path = {'human': human_protein_database_path,
                             'mouse': mouse_protein_database_path,
                             'rat': rat_protein_database_path,
                             'fish': fish_protein_database_path,
                             'fly': fly_protein_database_path,
                             'frog': frog_protein_database_path,
                             'mosquito': mosquito_protein_database_path}

    uniprot_mapping_database_path = 'database/AntibodyFilter_uniprot_id_mapping.csv'

    gene_information_database_path = 'database/Gene_Information.csv'

    fish_gene_id_mapping_database_path = 'database/GeneID_Mapping_Fish.csv'
    fly_gene_id_mapping_database_path = 'database/GeneID_Mapping_Fly.csv'
    frog_gene_id_mapping_database_path = 'database/GeneID_Mapping_Frog.csv'
    mosquito_gene_id_mapping_database_path = 'database/GeneID_Mapping_Mosquito.csv'
    worm_gene_id_mapping_database_path = 'database/GeneID_Mapping_Worm.csv'
    human_gene_id_mapping_database_path = 'database/GeneID_Mapping_Human.csv'
    mouse_gene_id_mapping_database_path = 'database/GeneID_Mapping_Mouse.csv'
    rat_gene_id_mapping_database_path = 'database/GeneID_Mapping_Rat.csv'

    gene_id_mapping_database_path = {'fish': fish_gene_id_mapping_database_path,
                                     'fly': fly_gene_id_mapping_database_path,
                                     'frog': frog_gene_id_mapping_database_path,
                                     'mosquito': mosquito_gene_id_mapping_database_path,
                                     'worm': worm_gene_id_mapping_database_path,
                                     'human': human_gene_id_mapping_database_path,
                                     'mouse': mouse_gene_id_mapping_database_path,
                                     'rat': rat_gene_id_mapping_database_path}

    human_ortholog_pair_best_database_path = 'database/Ortholog_Pair_Best_Human.csv'
    mouse_ortholog_pair_best_database_path = 'database/Ortholog_Pair_Best_Mouse.csv'
    rat_ortholog_pair_best_database_path = 'database/Ortholog_Pair_Best_Rat.csv'

    ortholog_pair_best_database_path = {'human': human_ortholog_pair_best_database_path,
                                        'mouse': mouse_ortholog_pair_best_database_path,
                                        'rat': rat_ortholog_pair_best_database_path}

    if mode == 'test':
        file_path = input_file_path
        output_path = f'test/{verified_target_species}_test_output_{start_time}.xlsx'
        err_output_path = f'test/{verified_target_species}_test_err_output_{start_time}.xlsx'
    else:
        file_path = input_file_path
        output_path = f'data/{verified_target_species}_output_{start_time}.xlsx'
        err_output_path = f'data/{verified_target_species}_err_output_{start_time}.xlsx'

    logger = Logger.Logger('CST_Phosphorylation_site_dataset_selected_March2024',
                           start_time,
                           'test' if mode == 'test' else 'prod',
                           f'get_antibody_for_{verified_target_species}')

    df = pd.read_excel(file_path)
    process_data(df,
                 output_path,
                 err_output_path,
                 logger,
                 verified_target_species,
                 uniprot_protein_database_path,
                 protein_database_path,
                 uniprot_mapping_database_path,
                 gene_information_database_path,
                 gene_id_mapping_database_path,
                 ortholog_pair_best_database_path)

    minutes, seconds = time_tracker.calculate_duration()

    logger.log(f'Total time taken: {minutes} minutes {seconds} seconds')


def process_data(df, output_path, err_output_path, logger, target_species,
                 uniprot_protein_database_path, protein_database_path,
                 uniprot_mapping_database_path, gene_information_database_path,
                 gene_id_mapping_database_path, ortholog_pair_best_database_path):
    process_count = 0
    updated_rows = []
    err_rows = []
    for _, row in df.iterrows():
        process_count += 1
        # Extracting specific columns by 1-based _: 1, 3, 5, 7, 10
        input_gene_name = row.iloc[0]  # 1st column
        input_acc_id = row.iloc[2]  # 3rd column
        site_info = row.iloc[4]  # 5th column
        input_species = row.iloc[6].lower()  # 7th column
        input_protein_slice = row.iloc[9].lower()  # 10th column
        cleaned_input_protein_slice = input_protein_slice.replace('_', '')
        specie_id_mapping = {
            'human': 9606,
            'mouse': 10090,
            'rat': 10116,
            'fly': 7227,
            'worm': 6239,
            'frog': 8364,
            'fish': 7955,
            'mosquito': 7165,
        }

        # Regular expression to extract amino acid and position from site_info
        match = re.match(r"([A-Za-z]+)(\d+)-p", site_info)
        if match:
            phospho_site_aa = match.group(1)
            input_target_position = int(match.group(2))
        else:
            logger.log(f'Error: Unable to extract amino acid and position from {site_info}')
            logger.log(f'Search failed')
            logger.blank()
            update_err(row,
                       err_rows,
                       'Unable to extract amino acid and position from site_info')
            continue

        logger.log(f'>>> [{process_count}] Current processing {input_species} gene: {input_gene_name}')

        logger.blank()
        logger.log(f'Starting Level 1 Search')
        logger.blank()

        # 蛋白质序列的获取
        # 一级搜索：通过uniprot_id获取
        input_protein_sequence = get_input_protein_sequence_by_uniprot_id(
            input_acc_id,
            cleaned_input_protein_slice,
            uniprot_protein_database_path,
            input_species)

        if input_protein_sequence:
            logger.log(f'Protein sequence found successfully for [{input_acc_id}] by uniprot_id')
            protein_source = 'uniprot_id'
            source_id = input_acc_id
        else:
            logger.log(f'Protein sequence not found for [{input_acc_id}] by uniprot_id')
            logger.log(f'Level 1 Search failed')
            logger.log(f'Starting Level 2 Search')
            # 二级搜索：通过RefSeq id获取
            input_protein_sequence = get_input_protein_sequence_by_refseq_id(
                input_acc_id,
                cleaned_input_protein_slice,
                protein_database_path,
                input_species)
            if input_protein_sequence:
                protein_source = 'refseq_id'
                source_id = input_acc_id
                logger.log(f'Protein sequence found successfully for [{input_acc_id}] by RefSeq id')
            else:
                logger.log(f'Protein sequence not found for [{input_acc_id}] by RefSeq id')
                logger.log(f'Level 2 Search failed')
                logger.log(f'Starting Level 3 Search')
                # 三级搜索：通过gene name获取
                input_protein_sequence, gene_id = get_input_protein_sequence_by_gene_id(
                    input_gene_name,
                    input_species,
                    cleaned_input_protein_slice,
                    gene_id_mapping_database_path,
                    protein_database_path)

                if input_protein_sequence and gene_id:
                    protein_source = 'gene_id'
                    source_id = gene_id
                    logger.log(f'Protein sequence found successfully for [{input_gene_name}] by {source_id}')
                else:
                    logger.log(f'Protein sequence not found for [{input_gene_name}] by gene name')
                    logger.log(f'Level 3 Search failed')
                    logger.blank()
                    update_err(row,
                               err_rows,
                               'Protein sequence not found from acc_id and gene name')
                    continue

        # 检查蛋白质序列中的目标氨基酸位置信息
        logger.blank()
        logger.log(f'Starting to check amino acid at position in protein sequence')
        logger.blank()

        # 检查 input_target_position 是否在有效范围内
        if 0 < input_target_position <= len(input_protein_sequence):
            # 比对目标氨基酸和指定位置的氨基酸
            if phospho_site_aa.lower() == input_protein_sequence[input_target_position - 1]:
                logger.log(f'Information Correct: Found amino acid "{phospho_site_aa}"'
                           f' at position "{input_target_position}" in protein sequence')
                correct_target_pos = input_target_position
            else:
                correct_target_pos = None
        else:
            correct_target_pos = None

        # 如果位置不正确或超出范围，尝试修复信息
        if correct_target_pos is None:
            logger.log(
                f'Information Incorrect: Amino acid "{phospho_site_aa}" not found '
                f'at position "{input_target_position}" in protein sequence')
            logger.log(
                f'Fixing information by searching for the amino acid "{phospho_site_aa}" in the protein sequence')

            correct_target_pos = find_actual_pos_by_slice(
                input_protein_sequence,
                input_protein_slice)

            if isinstance(correct_target_pos, str):
                logger.log(correct_target_pos)
                logger.log(f'Fixing failed')
                logger.blank()
                update_err(row, err_rows, 'Protein sequence found but not include target peptide sequence')
                continue

            logger.log(f'Corrected position: {correct_target_pos}')
        else:
            # 位置正确，无需修复
            pass

        logger.blank()
        logger.log(f'Starting to get gene id mapping')
        logger.blank()

        # 获取基因ID映射
        mapped_input_gene_id_list_by_uniprot = find_gene_id_mapping_by_uniprot_id(input_acc_id,
                                                                                  uniprot_mapping_database_path)

        if mapped_input_gene_id_list_by_uniprot:
            mapped_gene_id_counts = len(mapped_input_gene_id_list_by_uniprot)
            logger.log(f'{mapped_gene_id_counts} gene ids mapped successfully for {input_acc_id} by uniprot_id')
            mapped_input_gene_id_list = mapped_input_gene_id_list_by_uniprot
        else:
            logger.log(f'No gene id found for {input_acc_id} by uniprot_id')
            logger.log(f'Trying to find gene id by gene name')
            mapped_input_gene_id_list_by_gene_name = find_gene_id_mapping_by_gene_name(
                input_gene_name,
                input_species,
                gene_id_mapping_database_path)
            if mapped_input_gene_id_list_by_gene_name:
                mapped_gene_id_counts = len(mapped_input_gene_id_list_by_gene_name)
                logger.log(f'{mapped_gene_id_counts} gene ids mapped successfully for {input_gene_name} by gene name')
                mapped_input_gene_id_list = mapped_input_gene_id_list_by_gene_name
            else:
                logger.log(f'No gene id found for {input_gene_name} by gene name')
                logger.log(f'Search failed')
                logger.blank()
                update_err(row,
                           err_rows,
                           'Protein sequence found but no gene id found from uniprot_id and gene name')
                continue

        logger.blank()
        logger.log(f'Starting to search for orthologs in {target_species}')
        logger.blank()

        # 在目标物种中搜索同源基因
        target_species_id = specie_id_mapping.get(target_species)

        for idx_input_gene, gene_id in enumerate(mapped_input_gene_id_list):
            logger.log("---------------------------------------------")
            logger.log(f'[{process_count}-{idx_input_gene + 1}] Gene id: {gene_id}')
            logger.blank()
            logger.log(f'Searching for orthologs in {target_species}')
            logger.blank()
            # Search for orthologs in target species
            orthologs_target_species_gene_info_list = find_orthologs_by_gene_id(gene_id, input_species, target_species_id,
                                                                                ortholog_pair_best_database_path)
            if orthologs_target_species_gene_info_list:
                orthologs_count = len(orthologs_target_species_gene_info_list)
                logger.log(f'{orthologs_count} orthologs found in {target_species}')
            else:
                logger.log(f'No orthologs found in {target_species}')
                logger.log(f'Search failed')
                logger.blank()
                update_err(row,
                           err_rows,
                           f'Protein sequence and gene id found but no orthologs found in {target_species}')
                continue

            logger.blank()
            logger.log(f'Starting to get protein sequence by {target_species} gene id')
            logger.blank()

            for idx_target_species_gene, target_species_gene_info in enumerate(orthologs_target_species_gene_info_list):
                target_species_gene_id = target_species_gene_info[0]
                score = target_species_gene_info[1]
                best_score = target_species_gene_info[2]
                best_score_rev = target_species_gene_info[3]
                confidence = target_species_gene_info[4]

                logger.blank()
                logger.log(f'Searching for protein sequence by {target_species} gene id {target_species_gene_id}')
                logger.blank()

                target_species_specific_geneid, target_species_symbol = get_target_species_gene_id_info(
                    target_species_gene_id, gene_information_database_path)
                target_species_protein_sequence_list = get_target_species_protein_sequence_by_gene_id(
                    target_species_gene_id, protein_database_path, target_species)

                if target_species_protein_sequence_list:
                    if len(target_species_protein_sequence_list) > 1:
                        isoforms_count = len(target_species_protein_sequence_list)
                        logger.log(
                            f'{isoforms_count} protein isoforms found for {target_species} gene id: {target_species_gene_id}')
                        logger.blank()
                    else:
                        logger.log(
                            f'One protein sequence found successfully for {target_species} gene id: {target_species_gene_id}')
                else:
                    logger.log(f'Protein sequence not found for {target_species} gene id: {target_species_gene_id}')
                    logger.log(f'{target_species_gene_id} search failed')
                    logger.blank()
                    update_err(row,
                               err_rows,
                               f'Protein sequence, gene id and {target_species} orthologs found but no protein sequence found for {target_species} gene id')
                    continue

                for idx_protein, target_species_protein_info in enumerate(target_species_protein_sequence_list):
                    logger.log(".............................................")
                    logger.log(f'[{process_count}-{idx_target_species_gene + 1}-{idx_protein + 1}] '
                               f'Starting to align with input protein sequences')

                    output_protein_sequence = target_species_protein_info[0]
                    output_protein_source = target_species_protein_info[1]

                    aligned_input_seq, aligned_output_seq = align_protein_sequences(
                        input_protein_sequence,
                        output_protein_sequence)

                    aligned_target_common_pos = find_aligned_target_pos_in_input_sequence(
                        aligned_input_seq, correct_target_pos)

                    aligned_target_species_seq_slice = extract_slice(
                        aligned_output_seq, aligned_target_common_pos)

                    target_output_pos = find_actual_target_output_pos(
                        aligned_output_seq, aligned_target_common_pos)

                    aligned_input_seq_slice = extract_slice(aligned_input_seq, aligned_target_common_pos)

                    # logger.log(f'input_target_position: {input_target_position}')
                    # logger.log(f'aligned_target_common_pos: {aligned_target_common_pos}')
                    # logger.log(f'target_output_pos: {target_output_pos}')
                    # logger.blank()
                    # logger.log(f'input_protein_slice: {input_protein_slice}')
                    # logger.log(f'aligned_input_seq_slice: {aligned_input_seq_slice}')

                    logger.blank()
                    logger.log(f'{input_protein_slice}')
                    logger.log(f'{aligned_input_seq_slice}')
                    logger.log(f'{aligned_target_species_seq_slice}')
                    logger.blank()
                    # logger.log("==============================")
                    # logger.log(f'[{process_count}-{idx_target_species_gene + 1}-{idx_protein + 1}] Search Report:')
                    # logger.log("==============================")
                    # logger.blank()

                    (closest_p_pos,
                     count_p,
                     identity_15,
                     identity_11,
                     ) = get_report_info(
                        aligned_input_seq_slice,
                        target_output_pos,
                        aligned_target_species_seq_slice,
                        phospho_site_aa)

                    (best_window6,
                     highest_score6,
                     center_amino_match6,
                     side_scores_sum6,
                     offset6,

                     best_window7,
                     highest_score7,
                     center_amino_match7,
                     side_scores_sum7,
                     offset7,

                     best_window8,
                     highest_score8,
                     center_amino_match8,
                     side_scores_sum8,
                     offset8,

                     best_window9,
                     highest_score9,
                     center_amino_match9,
                     side_scores_sum9,
                     offset9,

                     best_window10,
                     highest_score10,
                     center_amino_match10,
                     side_scores_sum10,
                     offset10,

                     best_window11,
                     highest_score11,
                     center_amino_match11,
                     side_scores_sum11,
                     offset11
                     ) = get_window_info(
                        aligned_input_seq_slice,
                        aligned_target_species_seq_slice
                    )

                    # logger.log(f'Closest P position: {closest_p_pos}')
                    # logger.log(f'Count of P: {count_p}')
                    # logger.log(f'Identity for 15 length: {identity_15}')
                    # logger.log(f'Identity for 13 length: {identity_13}')
                    # logger.log(f'Identity for 11 length: {identity_11}')
                    # logger.blank()

                    # 数据整理
                    data = {
                        'order': process_count,
                        'protein_source': protein_source,
                        'source_id': source_id,
                        'target_phospho_pos': correct_target_pos,
                        'mapped_gene_id_from_uniprot_id': gene_id,
                        f'ortholog_{target_species}_gene_id': target_species_gene_id,
                        f'{target_species}_specific_geneid': target_species_specific_geneid,
                        f'{target_species}_symbol': target_species_symbol,
                        'DIOPT score': score,
                        'best_score': best_score,
                        'best_score_rev': best_score_rev,
                        'confidence': confidence,
                        f'{target_species}_protein_source': output_protein_source,
                        'aligned_input_seq_slice': aligned_input_seq_slice,
                        f'aligned_{target_species}_seq_slice': aligned_target_species_seq_slice,
                        f'target_{target_species}_pos': target_output_pos,
                        'closest_p_pos': closest_p_pos,
                        'count_p': count_p,
                        'identity for 15 length': identity_15,
                        'identity for 11 length': identity_11,

                        'best_window6': best_window6,
                        'highest_score6': highest_score6,
                        'target_phospho_site_match6': center_amino_match6,
                        'target_phospho_site_side_scores_sum6': side_scores_sum6,
                        'window6_offset': offset6,

                        'best_window7': best_window7,
                        'highest_score7': highest_score7,
                        'target_phospho_site_match7': center_amino_match7,
                        'target_phospho_site_side_scores_sum7': side_scores_sum7,
                        'window7_offset': offset7,

                        'best_window8': best_window8,
                        'highest_score8': highest_score8,
                        'target_phospho_site_match8': center_amino_match8,
                        'target_phospho_site_side_scores_sum8': side_scores_sum8,
                        'window8_offset': offset8,

                        'best_window9': best_window9,
                        'highest_score9': highest_score9,
                        'target_phospho_site_match9': center_amino_match9,
                        'target_phospho_site_side_scores_sum9': side_scores_sum9,
                        'window9_offset': offset9,

                        'best_window10': best_window10,
                        'highest_score10': highest_score10,
                        'target_phospho_site_match10': center_amino_match10,
                        'target_phospho_site_side_scores_sum10': side_scores_sum10,
                        'window10_offset': offset10,

                        'best_window11': best_window11,
                        'highest_score11': highest_score11,
                        'target_phospho_site_match11': center_amino_match11,
                        'target_phospho_site_side_scores_sum11': side_scores_sum11,
                        'window11_offset': offset11,
                    }

                    # Create a new dictionary for the row that includes both the original and new data
                    updated_row = row.to_dict()  # Convert the original row to a dictionary
                    updated_row.update(data)  # Update it with the new data

                    # Add the updated row to the result list
                    updated_rows.append(updated_row)

    # Create a new DataFrame with the enhanced rows
    update_df = pd.DataFrame(updated_rows)

    # Save the enhanced DataFrame to a new Excel file
    update_df.to_excel(output_path, index=False)

    err_df = pd.DataFrame(err_rows)
    err_df.to_excel(err_output_path, index=False)

    logger.log(f'> Finished processing all data')


def update_err(row, err_rows, err_info):
    err_row = row.to_dict()
    data = {
        'error': err_info
    }
    err_row.update(data)
    err_rows.append(err_row)


def get_input_protein_sequence_by_uniprot_id(uniprot_id, cleaned_input_protein_slice, uniprot_protein_database_path,
                                             input_species):
    # 使用 pandas 读取 CSV 文件
    df = pd.read_csv(uniprot_protein_database_path[input_species],
                     header=None, names=['uniprot_id', 'protein_sequence'])

    # 查找符合 uniprot_id 的行
    protein = df[df['uniprot_id'] == uniprot_id]

    if not protein.empty:
        # 获取 protein_sequence 并转换为小写
        protein_sequence = protein.iloc[0]['protein_sequence'].lower()
        # 检查是否包含 cleaned_input_protein_slice
        if cleaned_input_protein_slice in protein_sequence:
            return protein_sequence

    return None


def get_input_protein_sequence_by_refseq_id(refseq_id, cleaned_input_protein_slice,
                                            uniprot_protein_database_path, input_species):
    # 读取CSV文件, 没有标题行，因此指定header=None
    df = pd.read_csv(uniprot_protein_database_path[input_species], header=None,
                     names=['protein_acc', 'protein_sequence', 'geneid'])

    # 根据refseq_id查找对应的蛋白质序列
    protein_row = df[df['protein_acc'] == refseq_id]

    # 检查是否找到对应的蛋白质序列，并且该序列包含cleaned_input_protein_slice
    if not protein_row.empty and cleaned_input_protein_slice in protein_row.iloc[0]['protein_sequence']:
        return protein_row.iloc[0]['protein_sequence']
    else:
        return None


def get_input_protein_sequence_by_gene_id(gene_name, input_species, cleaned_input_protein_slice,
                                          gene_id_mapping_database_path, protein_database_path):
    print(input_species)
    gene_id_list = find_gene_id_mapping_by_gene_name(
        gene_name,
        input_species,
        gene_id_mapping_database_path)
    if gene_id_list:
        for gene_id in gene_id_list:
            df = pd.read_csv(protein_database_path[input_species], header=None,
                             names=['protein_version', 'protein_sequence', 'geneid'])
            # 根据geneid过滤数据
            protein = df[df['geneid'] == gene_id]
            if not protein.empty and protein['protein_sequence'].str.contains(cleaned_input_protein_slice).any():
                return protein['protein_sequence'].values[0], protein['geneid'].values[0]
    return None, None  # 保证在没有找到匹配时返回 (None, None)


def find_actual_pos_by_slice(protein_sequence, protein_slice):
    # Validate input length
    if len(protein_slice) < 8:
        return f"Error: '{protein_slice}' is less than 8 characters long."

    # Determine the actual index of the target in the cleaned version, accounting for underscores
    cleaned_index = 7 - protein_slice[:8].count('_')

    # Remove underscores from the slice for search purposes
    cleaned_protein_slice = protein_slice.replace('_', '')
    if not cleaned_protein_slice:
        return f"Error: '{protein_slice}' contains only underscores."

    # Attempt to find the cleaned slice in the protein sequence
    pos = protein_sequence.lower().find(cleaned_protein_slice.lower())
    if pos == -1:
        return f"'{protein_slice}' not found in '{protein_sequence}'."

    # Compute the actual position of the target character in the original sequence
    # Calculating the 1-based index position of the target character in the full protein sequence
    target_pos = pos + cleaned_index

    return target_pos + 1  # return 1-based index


def find_gene_id_mapping_by_uniprot_id(uniprot_id, uniprot_mapping_database_path):
    # 读取CSV文件
    df = pd.read_csv(uniprot_mapping_database_path, header=None, names=['uniprot_id', 'id_type', 'mapped_id'])

    # 筛选出uniprot_id匹配且id_type为'gene_id'的记录
    filtered_df = df[(df['uniprot_id'] == uniprot_id) & (df['id_type'] == 'gene_id')]

    # 提取mapped_id列表
    if not filtered_df.empty:
        gene_id_list = filtered_df['mapped_id'].tolist()
    else:
        gene_id_list = None

    return gene_id_list


def find_gene_id_mapping_by_gene_name(gene_name, input_species, gene_id_mapping_database_path):
    # 读取CSV文件，没有标题行，指定列名
    df = pd.read_csv(gene_id_mapping_database_path[input_species], header=None,
                     names=['geneid', 'idvalue2', 'species_specific_id'])

    # 过滤符合条件的记录
    filtered_records = df[(df['idvalue2'] == gene_name)]

    if not filtered_records.empty:
        gene_id_list = filtered_records['geneid'].tolist()
    else:
        gene_id_list = None

    return gene_id_list


def find_orthologs_by_gene_id(gene_id, input_species, target_tax_id, ortholog_pair_best_database_path):
    target_species_gene_info_list = []
    df = pd.read_csv(ortholog_pair_best_database_path[input_species], header=None,
                     names=['geneid1', 'species2', 'geneid2', 'score', 'best_score', 'best_score_rev', 'confidence'])

    # 将 gene_id 和 target_tax_id 以及 DataFrame 中的相关列都转换为字符串类型
    gene_id = str(gene_id)
    target_tax_id = str(target_tax_id)
    df['geneid1'] = df['geneid1'].astype(str)
    df['species2'] = df['species2'].astype(str)

    # 过滤满足条件的记录
    filtered_df = df[(df['geneid1'] == gene_id) &
                     (df['species2'] == target_tax_id) &
                     (df['confidence'].isin(['high', 'moderate']))]

    # 如果有匹配的记录，则将信息添加到列表中
    if not filtered_df.empty:
        for _, row in filtered_df.iterrows():
            target_species_gene_id = row['geneid2']
            score = row['score']
            best_score = row['best_score']
            best_score_rev = row['best_score_rev']
            confidence = row['confidence']
            target_species_gene_info_list.append(
                (target_species_gene_id, score, best_score, best_score_rev, confidence))
    else:
        target_species_gene_info_list = None

    return target_species_gene_info_list


def get_target_species_protein_sequence_by_gene_id(target_species_gene_id, protein_database_path, target_species):
    protein_sequence_list = []
    df = pd.read_csv(protein_database_path[target_species], header=None,
                     names=['protein_version', 'protein_sequence', 'geneid'])
    # 根据 gene_id 过滤数据
    filtered_df = df[df['geneid'] == target_species_gene_id]

    # 遍历过滤后的数据，提取蛋白质序列和版本信息
    for _, row in filtered_df.iterrows():
        seq = row['protein_sequence']
        source = row['protein_version']
        seq_info = (seq, source)
        protein_sequence_list.append(seq_info)

    return protein_sequence_list


def align_protein_sequences(input_protein_sequence, output_protein_sequence):
    # 使用pairwise2进行全局对齐
    alignments = pairwise2.align.globalms(
        input_protein_sequence.lower(),  # First sequence, converted to lowercase
        output_protein_sequence.lower(),  # Second sequence, converted to lowercase
        2,  # Match score
        -1,  # Mismatch penalty
        -0.5,  # Gap opening penalty
        -0.1  # Gap extension penalty
    )

    # 选择得分最高的对齐结果
    best_alignment = alignments[0]

    # 提取对齐后的两个序列
    aligned_seq1 = best_alignment[0]
    aligned_seq2 = best_alignment[1]

    return aligned_seq1, aligned_seq2


def find_aligned_target_pos_in_input_sequence(aligned_input_seq, position):
    aligned_position = 1
    actual_position = 0
    for i in range(len(aligned_input_seq)):
        if aligned_input_seq[i] != '-':
            actual_position += 1
        if actual_position == position:
            aligned_position = i + 1
            break

    return aligned_position


def extract_slice(sequence, pos):
    # pos是从1开始的位置，将其转换为从0开始的索引
    index = pos - 1

    # 计算片段的开始和结束位置
    start = max(index - 10, 0)
    end = min(index + 10 + 1, len(sequence))  # +1因为切片不包括结束位置

    # 提取片段
    segment = sequence[start:end]

    # 如果前面的字符不足7个，使用'-'补足
    if index - 10 < 0:
        segment = '-' * (10 - index) + segment

    # 如果后面的字符不足7个，使用'-'补足
    if index + 10 >= len(sequence):
        segment += '-' * (index + 10 - len(sequence) + 1)

    return segment


def find_actual_target_output_pos(aligned_output_seq, aligned_target_input_pos):
    cur_pos = 1
    actual_output_pos = 0
    for i in range(len(aligned_output_seq)):
        cur_pos += 1
        if aligned_output_seq[i] != '-':
            actual_output_pos += 1
        if cur_pos == aligned_target_input_pos:
            actual_output_pos += 1
            break

    return actual_output_pos


def get_report_info(aligned_input_seq_slice, target_output_phospho_pos, aligned_output_seq_slice, target_phospho_site):
    # 初始化相关变量
    identity_15 = 0
    identity_11 = 0
    count_p = 0
    closest_p_pos = None
    min_distance = float('inf')  # 用于比较找出最小距离
    phospho_site_options = {'y': ['y'], 's': ['s', 't'], 't': ['s', 't']}  # 可能的磷酸化位点
    valid_target_sites = phospho_site_options[target_phospho_site.lower()]

    # 中心位置索引（从0开始）
    center_index = len(aligned_output_seq_slice) // 2
    current_pos = target_output_phospho_pos

    # 辅助函数用于更新统计数据
    def update_stats(_, position_offset):
        nonlocal count_p, min_distance, closest_p_pos, current_pos
        current_pos += position_offset
        if aligned_output_seq_slice[i].lower() in valid_target_sites:
            count_p += 1
            distance = abs(current_pos - target_output_phospho_pos)
            if distance < min_distance:
                min_distance = distance
                closest_p_pos = current_pos

    # 检查和更新左侧序列
    for i in range(center_index - 1, -1, -1):
        if aligned_output_seq_slice[i] != '-':
            update_stats(i, -1)

    # 重置位置到中心
    current_pos = target_output_phospho_pos

    # 检查和更新右侧序列
    for i in range(center_index, len(aligned_output_seq_slice)):
        if aligned_output_seq_slice[i] != '-':
            update_stats(i, 1 if i > center_index else 0)

    # 计算identity_15, identity_11
    for i in range(len(aligned_output_seq_slice)):
        if aligned_output_seq_slice[i].lower() == aligned_input_seq_slice[i].lower() and aligned_output_seq_slice[
            i] != '-':
            if abs(i - center_index) <= 7:  # identity_15
                identity_15 += 1
            if abs(i - center_index) <= 5:  # identity_11
                identity_11 += 1

    return closest_p_pos, count_p, identity_15, identity_11


def get_window_info(aligned_input_seq_slice, aligned_output_seq_slice):
    best_window6, highest_score6, center_amino_match6, side_scores_sum6, offset6 = get_best_window(
        aligned_input_seq_slice.upper(),
        aligned_output_seq_slice.upper())

    best_window7, highest_score7, center_amino_match7, side_scores_sum7, offset7 = get_best_window(
        aligned_input_seq_slice.upper(),
        aligned_output_seq_slice.upper(),
        window_size=7)

    best_window8, highest_score8, center_amino_match8, side_scores_sum8, offset8 = get_best_window(
        aligned_input_seq_slice.upper(),
        aligned_output_seq_slice.upper(),
        window_size=8)

    best_window9, highest_score9, center_amino_match9, side_scores_sum9, offset9 = get_best_window(
        aligned_input_seq_slice.upper(),
        aligned_output_seq_slice.upper(),
        window_size=9)

    best_window10, highest_score10, center_amino_match10, side_scores_sum10, offset10 = get_best_window(
        aligned_input_seq_slice.upper(),
        aligned_output_seq_slice.upper(),
        window_size=10)

    best_window11, highest_score11, center_amino_match11, side_scores_sum11, offset11 = get_best_window(
        aligned_input_seq_slice.upper(),
        aligned_output_seq_slice.upper(),
        window_size=11
    )

    return (best_window6, highest_score6, center_amino_match6, side_scores_sum6, offset6,
            best_window7, highest_score7, center_amino_match7, side_scores_sum7, offset7,
            best_window8, highest_score8, center_amino_match8, side_scores_sum8, offset8,
            best_window9, highest_score9, center_amino_match9, side_scores_sum9, offset9,
            best_window10, highest_score10, center_amino_match10, side_scores_sum10, offset10,
            best_window11, highest_score11, center_amino_match11, side_scores_sum11, offset11)


def get_target_species_gene_id_info(gene_id, csv_file_path):
    # 读取CSV文件
    gene_data = pd.read_csv(csv_file_path, header=None, names=['geneid', 'symbol', 'species_specific_geneid'])

    # 查找对应的基因信息
    gene_info = gene_data[gene_data['geneid'] == gene_id].iloc[0]

    # 获取目标物种特定的geneid和symbol
    target_species_specific_geneid = gene_info['species_specific_geneid']
    target_species_symbol = gene_info['symbol']

    return target_species_specific_geneid, target_species_symbol


def get_best_window(aligned_input_seq_slice_length, aligned_output_seq_slice_length, window_size=6):
    # 根据中心点对称裁剪序列到window_size个氨基酸
    center_index = len(aligned_input_seq_slice_length) // 2
    aligned_input_seq_slice_length_window = aligned_input_seq_slice_length[
                                            center_index - window_size + 1:center_index + window_size]
    aligned_output_seq_slice_length_window = aligned_output_seq_slice_length[
                                             center_index - window_size + 1:center_index + window_size]

    # 载入BLOSUM62矩阵
    blosum62 = substitution_matrices.load("BLOSUM62")

    # 定义一个函数来计算匹配分数
    def score_amino_acids(aa1, aa2):
        if '-' in (aa1, aa2):  # 如果序列中有"-"字符，返回0分
            return 0.0
        else:
            score = blosum62[(aa1, aa2)]
            if aa1 == aa2:
                return 1.0
            elif score >= 0:
                return 0.8
            else:
                return 0.0

    # 定义一个函数来计算两个序列的匹配总分
    def calculate_total_score(seq1, seq2):
        total_score = 0
        for aa1, aa2 in zip(seq1, seq2):
            total_score += score_amino_acids(aa1, aa2)
        return total_score

    # 计算中心氨基酸左右两侧各一个氨基酸的匹配分数
    def calculate_side_scores(seq1, seq2):
        center = len(seq1) // 2
        left_score = score_amino_acids(seq1[center - 1], seq2[center - 1])
        right_score = score_amino_acids(seq1[center + 1], seq2[center + 1])
        return left_score + right_score

    # 定义一个函数找到窗口长度为指定值的最高分序列及其序列内容
    def find_highest_scoring_windows(seq1, seq2, windowSize):
        n = len(seq1)
        max_score = 0
        best_windows = []  # 存储所有最高得分的窗口和对应序列

        for i in range(n - windowSize + 1):
            window_score = calculate_total_score(seq1[i:i + windowSize], seq2[i:i + windowSize])
            if window_score > max_score:
                max_score = window_score
                best_windows = [(i, seq1[i:i + windowSize], seq2[i:i + windowSize])]
            elif window_score == max_score:
                best_windows.append((i, seq1[i:i + windowSize], seq2[i:i + windowSize]))

        return best_windows, max_score

    # 选择最对称的窗口，以序列的中心为基准
    def select_most_symmetric_window(best_windows, total_length, windowSize):
        sequence_center = total_length // 2
        min_diff = float('inf')
        most_symmetric_window = None
        final_offset = None

        for start, seq1, seq2 in best_windows:
            center_of_window = start + windowSize // 2  # 窗口中心索引
            current_diff = abs(center_of_window - sequence_center)
            if current_diff <= min_diff:
                min_diff = current_diff
                most_symmetric_window = (seq1, seq2)
                final_offset = center_of_window - sequence_center

        return most_symmetric_window, final_offset

    # 判断中心氨基酸是否相同
    if aligned_input_seq_slice_length_window[window_size - 1] == aligned_output_seq_slice_length_window[
        window_size - 1]:
        center_amino_match = "Yes"
    else:
        center_amino_match = "No"

    side_scores_sum = calculate_side_scores(
        aligned_input_seq_slice_length_window,
        aligned_output_seq_slice_length_window)

    highest_windows, highest_score = find_highest_scoring_windows(
        aligned_input_seq_slice_length_window,
        aligned_output_seq_slice_length_window,
        window_size)
    best_window, offset = select_most_symmetric_window(
        highest_windows,
        len(aligned_input_seq_slice_length_window),
        window_size)

    return best_window, highest_score, center_amino_match, side_scores_sum, offset


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python antibody_discovery.py [species](fly, frog, fish, mosquito, worm)")
        sys.exit(1)
    user_input_target_species = sys.argv[1]
    user_input_file_path = sys.argv[2] if len(sys.argv) > 2 else 'data/CST_Phosphorylation_site_dataset_selected_March2024.xlsx'
    try:
        main(user_input_target_species, user_input_file_path)
    except ValueError as e:
        print(e)
        sys.exit(1)
