# Protein structure analysis workflow


## Installation

The following commands install dependencies of the workflow in the _current_
directory within the `.pixi` folder. After installation, you cannot move the
folder without re-installling all the dependencies.

```bash
curl -fsSL https://pixi.sh/install.sh | bash
# ... cd <this repo>
pixi install
```

Because of the dependency on apptainer, this workflow only works on linux.

## Usage

Run `pixi run help` for the help page. Targets can be made using `pixi run make <target>`

## Example

In this example, we will compare the structures of yeast V-ATPase state 2 ([6O7W](https://www.rcsb.org/structure/6O7W)) to that of Citrus ([7UWB](https://www.rcsb.org/structure/7UWB))



In the terminal type

```bash
# you must type lowercase identifiers
pixi run make results/contact-map/6o7v results/contact-map/7uwb results/tm/7uw9-6o7w.score.txt results/plip/6o7v results/plip/7uwb--cores 5
```

To start the analysis with 5 paralel processes. This will download the structures and the programmes, and
run the analysis.

| ![](resources/contact-map-6o7v-vs-7uwb.png) |
|---------------------------------------------|
| Subunit interaction graph between `6o7v` (left) and `7uwb` (right |
