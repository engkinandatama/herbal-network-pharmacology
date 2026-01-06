"""
Command Line Interface
======================

Main CLI for the Network Pharmacology toolkit.
"""

import click
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn

from .config_loader import load_config

console = Console()


def print_header():
    """Print the toolkit header."""
    console.print("\n[bold green]═══════════════════════════════════════════════════════════════[/bold green]")
    console.print("[bold green]       Network Pharmacology Toolkit v1.0.0[/bold green]")
    console.print("[bold green]═══════════════════════════════════════════════════════════════[/bold green]\n")


@click.group()
@click.option(
    "--config", "-c",
    type=click.Path(exists=True),
    help="Path to project configuration YAML file"
)
@click.pass_context
def cli(ctx, config):
    """
    Network Pharmacology Toolkit - Analyze herbal medicine mechanisms.
    
    Use --config to specify a project configuration file.
    """
    ctx.ensure_object(dict)
    
    if config:
        ctx.obj["config"] = load_config(config)
    else:
        ctx.obj["config"] = None
    
    print_header()


@cli.command()
@click.pass_context
def info(ctx):
    """Display current project configuration information."""
    config = ctx.obj.get("config")
    
    if not config:
        console.print("[yellow]No config file specified. Use --config to load a project.[/yellow]")
        return
    
    table = Table(title="Project Configuration")
    table.add_column("Property", style="cyan")
    table.add_column("Value", style="green")
    
    table.add_row("Project Name", config.get("project.name", "N/A"))
    table.add_row("Plant", f"{config.plant_name} ({config.plant_latin_name})")
    table.add_row("Disease", config.disease_name)
    table.add_row("Known Compounds", str(len(config.known_compounds)))
    table.add_row("Data Directory", str(config.data_dir))
    
    console.print(table)


@cli.command()
@click.option("--source", "-s", 
              type=click.Choice(["all", "literature", "pubchem", "knapsack"]),
              default="all",
              help="Data source to collect from")
@click.pass_context
def collect(ctx, source):
    """Collect compound data from databases."""
    config = ctx.obj.get("config")
    
    if not config:
        console.print("[red]Error: No config file specified. Use --config option.[/red]")
        return
    
    console.print(f"[cyan]Collecting compounds for: {config.plant_name}[/cyan]")
    console.print(f"[cyan]Source: {source}[/cyan]\n")
    
    from .compounds.collector import CompoundCollector
    
    collector = CompoundCollector(config)
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        if source in ["all", "literature"]:
            task = progress.add_task("Loading compounds from literature...", total=None)
            compounds_lit = collector.collect_from_literature()
            progress.update(task, completed=True)
            console.print(f"  [green]✓[/green] Found {len(compounds_lit)} compounds from literature")
        
        if source in ["all", "pubchem"]:
            task = progress.add_task("Fetching from PubChem...", total=None)
            compounds_pubchem = collector.enrich_from_pubchem()
            progress.update(task, completed=True)
            console.print(f"  [green]✓[/green] Enriched with PubChem data")
        
        if source in ["all", "knapsack"]:
            task = progress.add_task("Searching KNApSAcK database...", total=None)
            compounds_knapsack = collector.search_knapsack()
            progress.update(task, completed=True)
            console.print(f"  [green]✓[/green] Found {len(compounds_knapsack)} compounds from KNApSAcK")
    
    # Save results
    output_file = collector.save_compounds()
    console.print(f"\n[green]Compounds saved to: {output_file}[/green]")


@cli.command()
@click.pass_context
def predict_targets(ctx):
    """Predict drug targets for compounds using SwissTargetPrediction."""
    config = ctx.obj.get("config")
    
    if not config:
        console.print("[red]Error: No config file specified. Use --config option.[/red]")
        return
    
    console.print(f"[cyan]Predicting targets for: {config.plant_name}[/cyan]\n")
    
    from .targets.predictor import TargetPredictor
    
    predictor = TargetPredictor(config)
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        task = progress.add_task("Loading compound data...", total=None)
        predictor.load_compounds()
        progress.update(task, completed=True)
        
        task = progress.add_task("Predicting targets via SwissTargetPrediction...", total=None)
        targets = predictor.predict_all()
        progress.update(task, completed=True)
        console.print(f"  [green]✓[/green] Predicted {len(targets)} target-compound pairs")
    
    output_file = predictor.save_targets()
    console.print(f"\n[green]Targets saved to: {output_file}[/green]")


@cli.command()
@click.pass_context
def get_disease_genes(ctx):
    """Retrieve disease-associated genes from databases."""
    config = ctx.obj.get("config")
    
    if not config:
        console.print("[red]Error: No config file specified. Use --config option.[/red]")
        return
    
    console.print(f"[cyan]Retrieving genes for: {config.disease_name}[/cyan]\n")
    
    from .disease.genes import DiseaseGeneCollector
    
    collector = DiseaseGeneCollector(config)
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        task = progress.add_task("Querying GeneCards...", total=None)
        genes_gc = collector.fetch_genecards()
        progress.update(task, completed=True)
        console.print(f"  [green]✓[/green] Found {len(genes_gc)} genes from GeneCards")
        
        task = progress.add_task("Querying DisGeNET...", total=None)
        genes_dg = collector.fetch_disgenet()
        progress.update(task, completed=True)
        console.print(f"  [green]✓[/green] Found {len(genes_dg)} genes from DisGeNET")
    
    # Merge and save
    all_genes = collector.merge_genes()
    output_file = collector.save_genes()
    console.print(f"\n[green]Disease genes saved to: {output_file}[/green]")
    console.print(f"[green]Total unique genes: {len(all_genes)}[/green]")


@cli.command()
@click.pass_context
def build_network(ctx):
    """Build and analyze PPI network."""
    config = ctx.obj.get("config")
    
    if not config:
        console.print("[red]Error: No config file specified. Use --config option.[/red]")
        return
    
    console.print("[cyan]Building protein-protein interaction network...[/cyan]\n")
    
    from .network.builder import NetworkBuilder
    from .network.analyzer import NetworkAnalyzer
    
    builder = NetworkBuilder(config)
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        task = progress.add_task("Finding common targets...", total=None)
        common = builder.find_common_targets()
        progress.update(task, completed=True)
        console.print(f"  [green]✓[/green] Found {len(common)} common targets")
        
        task = progress.add_task("Querying STRING database...", total=None)
        network = builder.build_ppi_network()
        progress.update(task, completed=True)
        console.print(f"  [green]✓[/green] Built network with {network.number_of_nodes()} nodes, {network.number_of_edges()} edges")
    
    # Analyze network
    analyzer = NetworkAnalyzer(config, network)
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        task = progress.add_task("Analyzing network topology...", total=None)
        hub_genes = analyzer.identify_hub_genes()
        progress.update(task, completed=True)
        console.print(f"  [green]✓[/green] Identified {len(hub_genes)} hub genes")
    
    # Save
    builder.save_network()
    analyzer.save_analysis()
    console.print(f"\n[green]Network analysis complete![/green]")


@cli.command()
@click.pass_context
def enrich(ctx):
    """Perform pathway enrichment analysis."""
    config = ctx.obj.get("config")
    
    if not config:
        console.print("[red]Error: No config file specified. Use --config option.[/red]")
        return
    
    console.print("[cyan]Performing enrichment analysis...[/cyan]\n")
    
    from .enrichment.analysis import EnrichmentAnalyzer
    
    analyzer = EnrichmentAnalyzer(config)
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        task = progress.add_task("Running KEGG pathway enrichment...", total=None)
        kegg = analyzer.kegg_enrichment()
        progress.update(task, completed=True)
        console.print(f"  [green]✓[/green] Found {len(kegg)} significant pathways")
        
        task = progress.add_task("Running GO enrichment...", total=None)
        go = analyzer.go_enrichment()
        progress.update(task, completed=True)
        console.print(f"  [green]✓[/green] Found {len(go)} GO terms")
    
    analyzer.save_results()
    console.print(f"\n[green]Enrichment analysis complete![/green]")


@cli.command()
@click.pass_context
def admet(ctx):
    """Predict ADMET properties for compounds."""
    config = ctx.obj.get("config")
    
    if not config:
        console.print("[red]Error: No config file specified. Use --config option.[/red]")
        return
    
    console.print("[cyan]Predicting ADMET properties...[/cyan]\n")
    
    from .admet.predictor import ADMETPredictor
    
    predictor = ADMETPredictor(config)
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        task = progress.add_task("Calculating drug-likeness...", total=None)
        predictor.predict_all()
        progress.update(task, completed=True)
    
    output_file = predictor.save_predictions()
    console.print(f"\n[green]ADMET predictions saved to: {output_file}[/green]")


@cli.command()
@click.option("--type", "-t",
              type=click.Choice(["network", "enrichment", "venn", "all"]),
              default="all",
              help="Type of visualization to generate")
@click.pass_context
def visualize(ctx, type):
    """Generate publication-ready figures."""
    config = ctx.obj.get("config")
    
    if not config:
        console.print("[red]Error: No config file specified. Use --config option.[/red]")
        return
    
    console.print(f"[cyan]Generating visualizations: {type}[/cyan]\n")
    
    from .visualization.plotter import Plotter
    
    plotter = Plotter(config)
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        if type in ["all", "network"]:
            task = progress.add_task("Plotting network...", total=None)
            plotter.plot_compound_target_network()
            plotter.plot_ppi_network()
            progress.update(task, completed=True)
            console.print("  [green]✓[/green] Network figures generated")
        
        if type in ["all", "venn"]:
            task = progress.add_task("Plotting Venn diagram...", total=None)
            plotter.plot_venn_diagram()
            progress.update(task, completed=True)
            console.print("  [green]✓[/green] Venn diagram generated")
        
        if type in ["all", "enrichment"]:
            task = progress.add_task("Plotting enrichment results...", total=None)
            plotter.plot_enrichment_barplot()
            plotter.plot_enrichment_bubble()
            progress.update(task, completed=True)
            console.print("  [green]✓[/green] Enrichment figures generated")
    
    console.print(f"\n[green]Figures saved to: {config.figures_dir}[/green]")


@cli.command()
@click.option("--skip-heavy", is_flag=True, help="Skip docking and MD (run in Colab)")
@click.pass_context
def run_all(ctx, skip_heavy):
    """Run the complete network pharmacology workflow."""
    config = ctx.obj.get("config")
    
    if not config:
        console.print("[red]Error: No config file specified. Use --config option.[/red]")
        return
    
    console.print("[bold green]Running complete workflow...[/bold green]\n")
    console.print(f"Plant: {config.plant_name}")
    console.print(f"Disease: {config.disease_name}\n")
    
    # Run all steps
    ctx.invoke(collect)
    ctx.invoke(predict_targets)
    ctx.invoke(get_disease_genes)
    ctx.invoke(build_network)
    ctx.invoke(enrich)
    ctx.invoke(admet)
    ctx.invoke(visualize)
    
    console.print("\n[bold green]═══════════════════════════════════════════════════════════════[/bold green]")
    console.print("[bold green]                   Workflow Complete!                          [/bold green]")
    console.print("[bold green]═══════════════════════════════════════════════════════════════[/bold green]")
    
    if skip_heavy:
        console.print("\n[yellow]Note: Docking and MD simulation skipped.[/yellow]")
        console.print("[yellow]Run notebooks/colab/ notebooks in Google Colab.[/yellow]")


def main():
    """Entry point for the CLI."""
    cli(obj={})


if __name__ == "__main__":
    main()
