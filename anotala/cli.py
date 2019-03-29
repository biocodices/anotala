from cement.core.foundation import CementApp
from cement.core.controller import CementBaseController, expose
from cement.utils.misc import init_defaults


defaults = init_defaults('anotala')

cache_argument = (
    ['--cache'],
    {'action': 'store', 'default': 'dict',
    'help': 'Cache option for annotator: mysql, redis, [dict]'}
)


class AppBaseController(CementBaseController):
    class Meta:
        label = 'base'


class ClinvarSearchController(CementBaseController):
    class Meta:
        label = 'clinvar'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            cache_argument,

            (['--phenos'],
             {'action': 'store', 'required': True,
              'help': ('Comma separated phenotypes. E.g: '\
                       '"Diabetes mellitus,Parkinson"')}),

            (['--clinsigs'],
             {'action': 'store',
              'default': 'pathogenic,likely pathogenic',
              'help': ('Comma separated clinical significances. Default: '\
                       '"pathogenic,likely pathogenic"')}),
        ]

    @expose(help='Search ClinVar variants for a list of phenotypes')
    def search(self):
        from anotala.recipes import search_clinvar_variants_for_phenotypes

        phenos = self.app.pargs.phenos.split(',')
        clinsigs = self.app.pargs.clinsigs.split(',')

        variants = search_clinvar_variants_for_phenotypes(
            pheno_terms=phenos,
            clinsigs=clinsigs,
            cache=self.app.pargs.cache,
        )

        self.app.render(variants)


class AnnotatePmidsController(CementBaseController):
    class Meta:
        label = 'pubmed'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            cache_argument,

            (['--pmids-file'],
             {'action': 'store', 'required': True,
              'help': 'File with the PMIDS to annotate, one per line.'})
        ]

    @expose(help='Anotate a list of PMIDs from a file')
    def annotate_pmids(self):
        from anotala import PubmedAnnotator

        fn = self.app.pargs.pmids_file
        annotator = PubmedAnnotator(cache=self.app.pargs.cache)

        with open(fn) as f:
            pmids = [line.strip() for line in f.readlines()]

        self.app.log.info('Read {} PMIDs from {}'.format(len(pmids), fn))

        annotations = annotator.annotate(pmids)

        self.app.render(annotations)


class AnotalaApp(CementApp):
    class Meta:
        label = 'anotala'
        config_defaults = defaults
        extensions = ['json', 'yaml']
        base_controller = 'base'
        handlers = [
            AppBaseController,
            ClinvarSearchController,
            AnnotatePmidsController,
        ]


with AnotalaApp() as app:
    app.run()
