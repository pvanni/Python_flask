from docx import Document


"""
Needed:
start_files directory
sequences directory
"""


report_file_name = "HOUSTON_report.docx"

doc = Document()
doc.add_heading('Report', 0)
para, pic = import_function(paired=True, sample_column='#SampleID')
report_writer(document=doc,
              heading="Importing data",
              paragraphs=para,
              pictures=None,)

doc.save(report_file_name)

para, pic, chosen_length = primer_removal(forward='CCTACGGGAGGCAGCAG', reverse='ATTACCGCGGCTGCTGG',
                                          threshold=0.90,
                                          no_trunc=True)

report_writer(document=doc,
              heading="Primer removal",
              paragraphs=para,
              pictures=pic)

doc.save(report_file_name)

para, pic, trunc_len = denoising(trunc_len=chosen_length, trim_left=0, paired=True, threshold=0.90)

report_writer(document=doc,
              heading="Denoising with DADA2",
              paragraphs=para,
              pictures=None)

doc.save(report_file_name)


para, pic = train_taxonomic_classifier(f_primer="CCTACGGGAGGCAGCAG",
                           r_primer="CCGTCAATTCMTTTRAGT",
                           trunc_len=trunc_len,
                           trim_left=0,
                           rep_file="/scratch/project_2001644/silva_138/seqs/silva138_representative_sequences.fasta",
                           tax_file="/scratch/project_2001644/silva_138/tax/silva138_taxonomy.tsv",
                           database_name="GreenGenes")

doc.save(report_file_name)

report_writer(document=doc,
              heading="Taxonomic classifier",
              paragraphs=para)

doc.save(report_file_name)

para, pic = filtering_table(f"denoising/{trunc_len}_dada2_output")

report_writer(document=doc,
              heading="Filtering the feature table",
              paragraphs=para)

doc.save(report_file_name)

para, pic = collapse_table("feature_table/feature_table.qza", 6)

report_writer(document=doc,
              heading="Collapsing",
              paragraphs=para)

doc.save(report_file_name)


para, pic = picrust2_table()

report_writer(document=doc,
              heading="PICRUSt2",
              paragraphs=para)

doc.save(report_file_name)


convert_to_ml(target_col='Delivery_mode',
              pos='C-section',
              neg='Vaginal')

