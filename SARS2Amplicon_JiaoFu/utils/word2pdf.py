import os
import win32com.client
from pathlib import Path
import click
import logging


logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)


def word2pdf(word, pdf):
    wdFormatPDF = 17
    inputFile = os.path.abspath(word)
    outputFile = os.path.abspath(pdf)
    word = win32com.client.Dispatch('Word.Application')
    doc = word.Documents.Open(inputFile)
    doc.SaveAs(outputFile, FileFormat=wdFormatPDF)
    doc.Close()
    word.Quit()

@click.command(context_settings={'help_option_names': ['-h','--help']})
@click.option('-w', '--workdir', required=True, type=click.Path(exists=True), help='WORD文件夹. 需要.docx后缀.')
def main(workdir):
    """Windows系统下运行, 批量WORD转PDF."""
    words = list(Path(workdir).resolve().glob('*.docx'))
    while words:
        word = words[0]
        pdf = str(word).replace('.docx','.pdf')
        try:
            word2pdf(word, pdf) # 有时候会报错
        except:
            logging.warning(f'WORD转化失败, 重试! - {word}')
        else:
            words.remove(word)

if __name__ == '__main__':
    main()
