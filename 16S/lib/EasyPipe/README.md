# README
微远基因所使用的流程构建模块
由于生产环境为单机，所以使用 parallel 进行并行计算

## 使用说明
### `Job`类
构建流程的基本原件，每个实例代表一个命令
```python
# 构建一个 Job 类的实例
job = Job("Prepare","dir_out")
job.set_workdir(f"{project.project_dir}/{prepare.name}")
template = "{perl} {script} -p {fastp} -l 150 -q 33 -r1 {rawdata}/{sample}_1.fq.gz -r2 {rawdata}/{sample}_1.fq.gz -o {dir_out}/{sample}"
for i in job.infos["samples"]:
    job.add_command(template.format(perl=project.softwares["perl"],
                                    script=os.path.join(project.bin,"fastp.pl"),
                                    fastp=project.softwares["fastp"],
                                    rawdata=project.infos["rawdata"],
                                    sample=i,
                                    dir_out=prepare.workdir))
```

### `ComplexJob`类
`Job` 类的复合体，每个实例代表一组命令，通过向其中添加 `Job` 类来生成复合命令

### `Pipe`类
`ComplexJob` 类的继承，每个实列代表一个流程，由基本的 `Job` 类 和 `ComplexJob` 类组成
```python

project = Pipe(project_name=project_name,
               project_dir=f"./{project_name}")

project.add(job)
```

