# Salve Mateus e Jess

O projetinho utiliza o renv para garantir a reprodutibilidade do ambiente R (semelhante ao uso de um requirements.txt no Python). Siga os passos abaixo para configurar tudo corretamente:

## Baixando o repositório localmente

`git clone https://github.com/C0rvito/DESeq2_Analysis.git`

## Entrando no projeto

Eu não sei porque, eu não entendi então não sei explicar...mas abram o projeto clicando no arquivo DESeq_Analysis.Rproj. Uma hora não entrei por aí e minha máquina quebrou

## Restaura as dependências (Caso não tenha (Imagino que você ja tenha Mateus))

```console
renv::restore()
```

## Rodando o pipeline

Modularizei o código, para termos uma reprodutibilidade melhor. Dá pra mudar o tipo de Câncer no main.R.
Rodar o main.R como se fosse python

# **ainda estou trabalhando nas utils\***
