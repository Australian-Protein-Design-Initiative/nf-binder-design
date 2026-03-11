# BoltzGen Nanobody Design

This example demonstrates how to design a nanobody binder for PDL1 using the BoltzGen workflow.

## Running the example

```bash
./run-local.sh
```

There are six nanobody scaffolds in the `nanobody_scaffolds` directory - five are from the [BoltzGen example](https://github.com/HannesStark/boltzgen/tree/main/example/nanobody_against_penguinpox), plus one additional scaffold (3eak, a humanized nanobody NbBCII10).

The `nanobody.yaml` configuration chooses scaffolds at random for each design.