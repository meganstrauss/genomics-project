# AMR Gene Atlas

Maps which antibiotic-resistance genes co-occur across bacterial genomes and
surfaces the cassettes driving multidrug resistance.

**Pipeline:** load gene catalogue + genomes → single-linkage gene clustering →
co-occurrence network weighted by lift (observed / expected co-occurrence).

```bash
python cli.py serve amr --port 8001   # from repo root
```

### API
| Method | Route                                 | Purpose                       |
| ------ | ------------------------------------- | ----------------------------- |
| POST   | `/api/pipeline`                       | ingest + cluster + network    |
| GET    | `/api/jobs/{id}`                      | poll the pipeline job         |
| GET    | `/api/network?min_lift=&min_count=`   | co-occurrence graph           |
| GET    | `/api/summary`                        | counts + carriage by class    |
