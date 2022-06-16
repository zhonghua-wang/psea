import typer
from psea import training
app = typer.Typer()


@app.command(help='prebuild require datasets, this step may take hours')
def build(nb_workers: int = 1):
    training.build(nb_workers=nb_workers)


@app.command(help='Download or update the HPO database')
def download_hpo(override: bool = False):
    from psea import hpo_config
    conf = hpo_config.HPOConfig()
    conf.init_hpo(override=override)


@app.command(help='Priorities candidates Gene for given HPO terms')
def predict(hpo_terms: str, top: int = 10, output_file: str = None):
    pass


if __name__ == '__main__':
    app()
