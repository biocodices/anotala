from cement.core.foundation import CementApp
from cement.core.controller import CementBaseController, expose
from cement.utils.misc import init_defaults

from anotamela.recipes import search_clinvar_variants_for_phenotypes


defaults = init_defaults('anotamela')


class BaseController(CementBaseController):
    class Meta:
        label = 'base'
        arguments = [
            (['--phenos'],
             {'action': 'store', 'required': True,
              'help': ('Comma separated phenotypes. E.g: '\
                       '"Diabetes mellitus,Parkinson"')}),

            (['--clinsigs'],
              {'action': 'store', 'required': True,
               'default': 'pathogenic,likely pathogenic',
               'help': ('Comma separated clinical significances. Default: '\
                        '"pathogenic,likely pathogenic"')}),

            (['--cache'],
             {'action': 'store', 'default': 'dict',
              'help': 'Cache option for annotator: mysql, redis, [dict]'}),
        ]

    @expose(help='Search ClinVar variants for a list of phenotypes')
    def search_clinvar(self):
        phenos = self.app.pargs.phenos.split(',')
        clinsigs = self.app.pargs.clinsigs.split(',')

        variants = search_clinvar_variants_for_phenotypes(
            pheno_terms=phenos,
            clinsigs=clinsigs,
            cache=self.app.pargs.cache,
        )

        app.render(variants)


class AnotamelaApp(CementApp):
    class Meta:
        label = 'anotamela'
        config_defaults = defaults
        extensions = ['json', 'yaml']
        base_controller = 'base'
        handlers = [BaseController]


with AnotamelaApp() as app:
    app.run()
