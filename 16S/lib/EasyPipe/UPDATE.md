# UPDATE LOG 

## v4.0.0
专为微远公卫生信开发的单机版本，此版本后续需同之前的 pbs 版本统一接口
### fix
- None
### update
- 重写 Job、ComplexJob、Pipe 类，使之适配单机
- 移除了 utils.py 中的 shell_run

## v3.1.0
重新整理了一些内部逻辑，现在构建流程最好使用对Pipe对象的继承去做
### fix
- None
### update
- None

## v3.0.0
代码重构，现在使用堆叠的方式去构建流程
### fix
- None
### update
- None


## beta 2.0.0
### fix
- None
### update
- Add some doc for module description
- Add UPDATE.md
- Add timeinfo decorator