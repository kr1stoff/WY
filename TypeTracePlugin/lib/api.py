#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/12/19 15:48
# @Last Modified by:   Ming
# @Last Modified time: 2023/12/19 15:48
import json
import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import urljoin

import requests

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def downloader(url, f_out, timeout=300, mode='w'):
    """
    下载内容到指定的文件

    :param url: The download url
    :param f_out: The output file name
    :param timeout: The timeout for download
    :param mode: The write mode
    """
    try:
        response = requests.get(url, timeout=timeout, stream=True)
        if response.status_code == 200:
            if str(f_out).endswith('.gz'):
                OUT = open(f_out, 'wb')
            else:
                OUT = open(f_out, mode)
            for chunk in response.iter_content(chunk_size=128):
                OUT.write(chunk)
            OUT.close()
            return 0
        else:
            return (url, f_out)
    except Exception as e:
        logger.debug(e)
        return (url, f_out)


def multi_thread_downloader(info: dict[str, Path], f_fail: Path, threads=5, force=True, mode='w'):
    """
    多线程下载

    :param info: A dict with url as key and file path as value
    :param f_fail: A file to store the fail info
    :param threads: The threads to use
    :param force: Whether download the file that already exists
    :param mode: The write mode
    """
    tasks = []
    fail_tasks = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for f_url, f_out in info.items():
            if not f_out.exists() or force:
                tasks.append(executor.submit(downloader, f_url, f_out, mode=mode))
        for future in as_completed(tasks):
            exec_res = future.result()
            if exec_res != 0:
                fail_tasks.append(future.result())

    if len(fail_tasks) > 0:
        logger.warning(f"Write the fail content to {f_fail}")
        with open(f_fail, 'w') as OUT:
            for fail_task in fail_tasks:
                print(*fail_task, sep='\t', file=OUT)
    else:
        f_fail.unlink(missing_ok=True)


class EnteroBase(object):
    """
    EnteroBase数据库下载工具
    """

    def __init__(self):
        """
        Init the object
        """
        self.baseurl = "https://enterobase.warwick.ac.uk/schemes/"

    def list(self):
        """
        List the available schemes in the EnteroBase database
        """
        session = requests.get(self.baseurl)
        self.schemes = []
        for i in re.findall(r'<a href=\"(\S+)/\">', session.text):
            if i != "..":
                self.schemes.append(i)
        return self.schemes

    def downloader(self, url, f_out, timeout=300):
        """
        下载内容到指定的文件

        :param url: The download url
        :param f_out: The output file name
        :param timeout: The timeout for download
        """
        try:
            response = requests.get(url, timeout=timeout, stream=True)
            if response.status_code == 200:
                if str(f_out).endswith('.gz'):
                    OUT = open(f_out, 'wb')
                else:
                    OUT = open(f_out, 'w')
                for chunk in response.iter_content(chunk_size=128):
                    OUT.write(chunk)
                OUT.close()
                return 0
            else:
                return (url, f_out)
        except Exception as e:
            return (url, f_out)

    def download(self, scheme, d_out, threads=5, force=False):
        """
        Download the available schemes to the target directioy

        :param scheme: The scheme you want to download
        :param d_out: The output dir
        :param threads: The thread number to use
        :param force: Whether redownload the files exists
        """
        self.list()
        if scheme not in set(self.schemes):
            raise ValueError(scheme)
        else:
            url = urljoin(self.baseurl, scheme)
            file_to_download = []
            for i in re.findall(r'<a href=\"(\S+)\">', requests.get(url).text):
                if i != "../":
                    file_to_download.append(i)

        tasks = []
        fail_tasks = []
        with ThreadPoolExecutor(max_workers=threads) as executor:
            for i in file_to_download:
                f_url = urljoin(self.baseurl, f"{scheme}/{i}")
                f_out = d_out.joinpath(i)
                if not f_out.exists() or force:
                    tasks.append(executor.submit(self.downloader, f_url, f_out))
            for future in as_completed(tasks):
                exec_res = future.result()
                if exec_res != 0:
                    fail_tasks.append(future.result())

        if len(fail_tasks) > 0:
            with open(d_out.joinpath("fail.txt"), 'w') as OUT:
                for fail_task in fail_tasks:
                    print(*fail_task, sep='\t', file=OUT)

    def download_fail(self, f_fail: Path, threads=5, force=True):
        """
        Download the file in the fail.txt

        :param f_fail: The file of fail.txt
        :param threads: The threads to download content
        :param force: If true, overwrite existing file
        """
        info_download = {}
        with open(f_fail, 'r') as IN:
            for line in IN:
                arr = line.strip().split("\t")
                info_download[arr[0]] = Path(arr[1])
        multi_thread_downloader(info_download, f_fail, threads, force, mode='wb')


class Pasteur(object):
    """
    Pasteur网站数据库工具
    """

    def __init__(self):
        """
        Init the object
        """
        self.base_url = "https://bigsdb.pasteur.fr/api/"

    def list_scheme(self):
        """
        下载并展示数据库内包含的物种及Scheme
        """
        res = {}
        response = requests.get(self.base_url, timeout=300)
        if response.ok:
            resources = response.json()
        else:
            return response.status_code

        for resource in resources:
            for database in resource["databases"]:
                if "seqdef" in database["href"]:
                    name = database["name"]
                    logger.info(f"{name}")
                    res[name] = {}
                    definitions_page = requests.get(database["href"], timeout=300)
                    if definitions_page.ok:
                        schemes_url = definitions_page.json()["schemes"]
                        info_schemes = requests.get(schemes_url, timeout=300).json()
                        for scheme in info_schemes["schemes"]:
                            desc = scheme["description"]
                            res[name][desc] = scheme["scheme"]
                    else:
                        logger.error(f"{name}: {definitions_page.status_code}")
        return res

    def download_scheme(self, database, scheme_id, d_out: Path, prefix, threads, force):
        """
        下载指定的Scheme及等位基因序列

        :param database: The database name such as pubmlst_arcobacter_seqdef
        :param scheme_id: The scheme id such as 1
        :param d_out: The output dir
        :param prefix: The prefix for scheme out put
        :param threads: The max threads to use when download the allele sequence
        :param force: Whether download the exists file again
        """
        d_out = d_out.absolute()
        d_out.mkdir(exist_ok=True, parents=True)

        # 下载info信息页
        info_url = urljoin(self.base_url, f"db/{database}/schemes/{scheme_id}")
        f_info = d_out.joinpath(f"info.json")
        if not f_info.exists() or force:
            logger.info(f"Download the scheme info to {f_info}")
            info = requests.get(info_url, timeout=300).json()
            with open(f_info, 'w') as OUT:
                json.dump(info, OUT, indent=2)
        else:
            with open(f_info, 'r') as IN:
                info = json.load(IN)

        # 下载scheme
        scheme_url = urljoin(self.base_url, f"db/{database}/schemes/{scheme_id}/profiles_csv")
        f_scheme = d_out.joinpath(f"{prefix}.tsv")
        if not f_scheme.exists() or force:
            logger.info(f"Download the scheme to {f_scheme}")
            res = downloader(scheme_url, f_scheme)
            if res != 0:
                logger.warning(f"No scheme find, Pleast check the url {scheme_url}")

        # 下载等位基因序列
        logger.info(f"Download the allele sequence to {d_out}")
        num_loci = len(info["loci"])
        threads = min(num_loci, threads)
        f_fail = d_out.joinpath("fail.txt")
        info_download = {}
        for loci in info["loci"]:
            name = loci.strip().split('/')[-1]
            url_allele = f"{loci}/alleles_fasta"
            f_allele = d_out.joinpath(f"{name}.fasta")
            info_download[url_allele] = f_allele
        multi_thread_downloader(info_download, f_fail, threads, force, mode='wb')

    def download_fail(self, f_fail: Path, threads=5, force=True):
        """
        Download the file in the fail.txt

        :param f_fail: The file of fail.txt
        :param threads: The threads to use
        :param force: Whether download the exists file again
        """
        info_download = {}
        with open(f_fail, 'r') as IN:
            for line in IN:
                arr = line.strip().split("\t")
                info_download[arr[0]] = Path(arr[1])
        multi_thread_downloader(info_download, f_fail, threads, force, mode='wb')
