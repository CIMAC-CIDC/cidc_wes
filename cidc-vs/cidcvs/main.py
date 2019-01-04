""" Visualization of WES harmonization 2
Usage: cidc-vs /path/to/files /path/to/output

Input:

    Path of folder for all files

    Path of folder for output figures

Files nomenclature:

    Center_source_sampleID_tools_feature
    Options:
        Center: Broad, MDACC
        tools: Sentieon, GISTIC, Pyclone, Facets
        source: FFPE, Frozen
        feature: Mutation.Germline, Mutation.Somatic, Clonality, Purity, CNV

"""

from .preprocess import *
from .definition import *
from .plot import *
import warnings
from matplotlib.backends.backend_pdf import PdfPages

BASE_DIR = os.path.dirname(os.path.realpath(__file__))+'/'
REF = os.path.join(BASE_DIR, 'data/REF')
class Compare(object):
    __slots__ = [
                "fmt","DPI","purity_dir", "clonality_dir", "cnv_dir", "out_dir",
                "somatic_mutation_maf_dir", "somatic_mutation_vcf_dir",
                "indel_mutation_maf_dir", "sv_mutation_maf_dir", "ref_fasta", "cancer",
                ]

    def __init__(self, purity_dir,clonality_dir,cnv_dir,out_dir,
                somatic_mutation_maf_dir,somatic_mutation_vcf_dir,
                 indel_mutation_maf_dir, sv_mutation_maf_dir, ref_fasta,cancer):
        self.cancer = cancer
        self.ref_fasta = ref_fasta
        self.purity_dir = purity_dir
        self.clonality_dir = clonality_dir
        self.cnv_dir = cnv_dir
        self.out_dir = out_dir
        self.somatic_mutation_maf_dir = somatic_mutation_maf_dir
        self.somatic_mutation_vcf_dir = somatic_mutation_vcf_dir
        self.indel_mutation_maf_dir = indel_mutation_maf_dir
        self.sv_mutation_maf_dir = sv_mutation_maf_dir
        self.DPI = 100
        self.fmt='png'
   
    # purity
    def purity(self):
        facet_result = concatPurity(base_path=self.purity_dir)
        purity = pd.pivot_table(facet_result, index=[
                                'ID', 'Source'], columns='Center', values='purity').reset_index()
        ploidy = pd.pivot_table(facet_result, index=[
                                'ID', 'Source'], columns='Center', values='ploidy').reset_index()
        plt.figure(figsize=(6, 12))
        ax = plt.subplot(211)

        qqPlot(df=purity ,x='Broad', y='MDA', factor='Source', ax=ax, label='ID')
        ax.set_title('Estimated Tumor Purity')
        ax = plt.subplot(212)
        
        qqPlot(df=ploidy, x='Broad', y='MDA',  factor='Source', ax=ax, label='ID')
        ax.set_title('Estimated Tumor Ploidy')

        plt.subplots_adjust(hspace=0.3)
        plt.savefig('{0}/purity.{1}'.format(self.out_dir,self.fmt), dpi=self.DPI)
        plt.clf()

    # clonality
    def clonality(self):
        pyclone_out = concatClonality(base_path=self.clonality_dir)
        plt.figure(figsize=(6, 10))
        ax = plt.subplot(211)
        barPlot(x="ID", y="Clonality", hue="Center", 
                df=pyclone_out.loc[pyclone_out.Source == 'FFPE', :],
                title='Clonality of FFPE Samples ',
                ax=ax, hue_order=['MDA', 'Broad'])
        plt.legend(loc=(.8, -0.3))

        ax = plt.subplot(212)
        barPlot(x="ID", y="Clonality", hue="Center", 
                df=pyclone_out.loc[pyclone_out.Source == 'Frozen', :],
                title='Clonality of Frozen Samples ',
                ax=ax, hue_order=['MDA', 'Broad'])
        ax.legend_.remove()
        plt.subplots_adjust(top=0.8, hspace=0.5)
        plt.savefig('{0}/clonality.{1}'.format(self.out_dir,self.fmt), dpi=self.DPI)
        plt.clf()


    # CNV 
    def cnv(self):
        cnv_cor,metadata = concatCNV(base_path=self.cnv_dir)
        colors = {'Center': ['peachpuff', 'lightblue'],
                'Source': ['aliceblue', 'lightgray']}
        heatMap(df=cnv_cor,col_colors=metadata,col_color_dict=colors,figsize=(12,12))
        plt.subplots_adjust(right=0.8)
        plt.savefig('{0}/cnv.{1}'.format(self.out_dir,self.fmt), dpi=self.DPI)
        plt.clf()

    # somatic mutation
    def somatic_mutation(self):
        somatic_maf,label = calJaccardMtrx(base_path=self.somatic_mutation_maf_dir)
        label = label.map(dict(zip( 
                                sorted(label.unique()),
                                ['peachpuff', 'lightblue']
                                )))

        plt.figure(figsize=(7, 5))
        ax = plt.subplot(111)
        dendrogramPlot(df= 1 - somatic_maf,label=label,ax=ax)
        plt.subplots_adjust(right=.7, bottom=0.2)
        plt.savefig(
            '{0}/somatic_mutation_cluster_tnsnv.{1}'.format(self.out_dir,self.fmt), dpi=self.DPI)
        plt.clf()

        ## Mutation profile
        with PdfPages('{}/somatic_mutation_profile_tnsnv.pdf'.format(self.out_dir)) as pdf:
            for sampleID,tri_mtrx in iterMaf(base_path=self.somatic_mutation_maf_dir, ref_fasta=self.ref_fasta):

                ref_mtrx = pd.concat([pd.read_table("{0}/{1}.mtrx".format(REF,x)) for x in self.cancer],axis=0)
                tri_mtrx['Group'] = sampleID
                fig = plot96Mtrx(df=pd.concat(
                    [tri_mtrx, ref_mtrx]), rows='Group')
                pdf.savefig(fig,dpi=self.DPI)
                plt.clf()
      
            

        
    def indel_mutation(self):
        indel_maf, label = calJaccardMtrx(
            base_path=self.indel_mutation_maf_dir,mut_type='indel')
        
        label = label.map(dict(zip(
            sorted(label.unique()),
            ['peachpuff', 'lightblue']
        )))

        plt.figure(figsize=(7, 5))
        ax = plt.subplot(111)
        dendrogramPlot(df=1 - indel_maf, label=label, ax=ax)
        plt.subplots_adjust(right=.7, bottom=0.2)
        plt.savefig(
            '{0}/indel_mutation_cluster.{1}'.format(self.out_dir, self.fmt), dpi=self.DPI)
        plt.clf()

    def sv_mutation(self):
        # TODO Waiting for result running on myc server, and then convert vcf to maf
        # Store in somatic_mutation/sv_maf

        sv_maf, label = calJaccardMtrx(
            base_path=self.sv_mutation_maf_dir, mut_type='sv')
        label = label.map(dict(zip(
            sorted(label.unique()),
            ['peachpuff', 'lightblue']
        )))

        plt.figure(figsize=(7, 5))
        ax = plt.subplot(111)
        dendrogramPlot(df=1 - sv_maf, label=label, ax=ax)
        plt.subplots_adjust(right=.7, bottom=0.2)
        plt.savefig(
            '{0}/sv_mutation_cluster.{1}'.format(self.out_dir, self.fmt), dpi=self.DPI)
        plt.clf()
        

   


    
    
    
