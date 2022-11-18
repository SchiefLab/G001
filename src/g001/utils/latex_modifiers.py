# -*- coding: utf-8 -*-
from pathlib import Path
import re
import shutil
import subprocess
from typing import Tuple, List

import docx2pdf
from PyPDF2 import PdfFileMerger


def latex2pdf(path: Path, outpath: Path, caption_prefix: str = None, cleanheaders: bool = True) -> None:
    """
    Takes a complete document latex files and will strip it of it's pages
    numbers and it's figure/table captions prefix in place of your own. The
    DOCX files is then converted into a PDF
    Parameters
    ----------
    path : Path
        path to the latex file
    outpath : Path
        path to resulting pdf file
    caption_prefix : str
        prefix you want after "Table" i.e. "Table S1" if "S1" is the input
    """
    path_fixed = path.with_stem(path.stem + "-tmp")
    # Add correct caption and remove any pages numbers
    with open(path) as infile:
        tmp = (
            infile.read()
            .replace(r"\end{document}", "\\thispagestyle{empty}\n\\pagenumbering{gobble}\n\\end{document}")
            .replace("\\begin{table}", "\\captionsetup{labelformat=empty}\n\\begin{table}")
            .replace("\\begin{document}", "\\usepackage{caption}\n\\begin{document}")
        )
        final = []
        for line in tmp.split("\n"):
            # Uses colon as an anchor for caption!
            if "\\caption[" in line and caption_prefix and ":" in line:
                # Real caption replacement -- if inner caption exists that is the one actually used by LaTex
                for outter_caption in re.finditer(r"caption\[[\(\w\s\)]*\:", line):
                    start, end = outter_caption.span()
                    line = line[:start] + f"caption[{caption_prefix}" + line[end:]
                # Inner priority caption
                for inner_caption in re.finditer(r"\}[\((\w\s\)]*\:", line):
                    start, end = inner_caption.span()
                    if cleanheaders:
                        caption_prefix = "".join(string_profiler(caption_prefix))
                    line = line[:start] + "}" + f"{caption_prefix}" + line[end:]
            final.append(line)
        final = "\n".join(final)
        with open(path_fixed, "w") as f:
            f.write(final)
    # we add igore errors "nonstopmode" because it will fail for shitty reasons.
    _ = subprocess.run(
        ["xelatex", "-interaction=nonstopmode", "-output-directory", path.parent, path_fixed],
        capture_output=True,
    )
    # rename tmp pdf to real pdf name
    path_fixed.with_suffix(".pdf").rename(path.with_suffix(".pdf"))
    # del all the junk files
    for file in path.parent.glob("*-tmp*"):
        file.unlink()


def make_tables(
    srcpath: Path,
    outpath: Path,
    pdf_order: list[Tuple[str, str, str]],
    cleanheaders: bool = False,
    paper_numnber: int = 1,
) -> None:
    """
    Takes a list of tables and creates a pdf with them in the order you want
    Parameters
    ----------
    outpath : Path
        path to resulting pdf file
    pdf_order : Tuple[str, str]
        ('S1', '/path/to/S1.tex')
    Raises
    ------
    FileNotFoundError
        If the file is not found
    ValueError
        If the file is not a valid file type, must be (pdf, tex, docx)
    """
    merger = PdfFileMerger()
    for caption_prefix, report_table, filename in pdf_order:
        path = srcpath / filename
        if not path.exists():
            raise FileNotFoundError(f"{path} does not exist")
        if path.suffix == ".docx":
            # TODO: we cant update the docx file so any modifcations need to be done manually
            docx2pdf.convert(str(path), str(path.with_suffix(".pdf")), keep_active=True)
        elif path.suffix == ".tex":
            report_table = f" ({report_table})" if report_table else ""
            prefix = f"Table {caption_prefix}{report_table}:"  # because im lazy and dont want to type this
            latex2pdf(path, path.with_suffix(".pdf"), prefix, cleanheaders)
        elif path.suffix != ".pdf":
            raise ValueError(f"{path} is not a valid file type")
        pdf_path = path.with_suffix(".pdf")
        merger.append(str(pdf_path))  # old pdf merger cannot handle Path objects
        shutil.copy2(pdf_path, outpath / f"pdfs/{caption_prefix + '.pdf'}")

    merger.write(str(srcpath / "all-sup-tables.pdf"))
    merger.close()
    process = subprocess.run(
        ["xelatex", "-output-directory", outpath, srcpath / f"SupplementalTablesStagingPaper.tex"],
        capture_output=False,
    )
    process.check_returncode()
    (srcpath / "all-sup-tables.pdf").unlink()  # A confusing file to keep so just delete it.


def string_profiler(
    string: str,
    start_delimiter: str = "(",
    end_delimiter: str = ")",
    remove: bool = True,
    keep_delimiter: bool = True,
) -> List[str]:
    """
    Seperates strings fragements into list based on the start and end delimiters
    Args:
        string: complete string you want to be broken up based on start and stop delimiters given
        start_delimiter: delimiter element to start
        end_delimiter: delimiter elemtent to end
        remove: decide whether or not to keep strings inside the delimiters
    Returns:
        List[str]: list of strings that are split at start and end delimiters given and whether
            or not you want to remove the string inside the delimiters
    Tests:
        long = '(life is is good) love world "(blah) blah" "here I am" once again "yes" blah '
        print(string_profiler(long))
        null = ''
        print(string_profiler(null))
        short = '(life love) yes(and much more)'
        print(string_profiler(short))
        short = 'yes "life love"'
        print(string_profiler(short))
    """
    outer_index = 0  # stepper for outer delimier string elements
    inner_index = 0  # stepper for inner delimier string elements
    curr_index = 0  # actual index of the current element in the string
    string_list = []  # string broken up into individual elements whenever a start and end delimiter is hit
    outer_string = ""  # string tracked while outside the delimiters
    inner_string = ""  # string tracked while inside the delimiters

    for outer_index in range(len(string)):
        # Actual pointer position (inner delimiter counter + outer delimiter counter)
        curr_index = inner_index + outer_index
        # Close once acutal index is at the end
        # NOTE: outer_index will keep going till end regardless of hitting a delimiter and adding to inner stepper.
        if curr_index == len(string):
            break
        ### DELIMITER HIT ###
        if string[curr_index] == start_delimiter:
            # If we his a delimiter, collect the string previous to that as an element; flush
            if outer_string:
                # Option: .extend(outer_string.strip().split()) | If you want every word seperate. Maybe an option?
                string_list.append(outer_string.strip())
                outer_string = ""
            for j in range(curr_index + 1, len(string)):
                # Stepper that is pushed while in inner delimiter string.
                inner_index += 1
                # Once we his the end delimiter, stop iterating through the inner delimiter string
                if string[j] == end_delimiter:
                    break
                # String inside delimiters
                inner_string += string[j]
            # If you want the string inside the delimiters
            if not remove:
                if keep_delimiter:
                    inner_string = start_delimiter + inner_string + end_delimiter
                string_list.append(inner_string)
            # inner delimiter string restart
            inner_string = ""
        # String outside of the delimiters
        else:
            outer_string += string[curr_index]
        # End delimiter is either nested or not the real target; should ignore
        if string[curr_index] == end_delimiter:
            if string_list and outer_string:
                string_list[-1] += outer_string
                outer_string = ""
    # In case of not hiting a delimiter at the end of the string, collect the remaining outer delimiter string
    # Option: .extend(outer_string.strip().split()) | If you want every word seperate. Maybe an option?
    if outer_string:
        string_list.append(outer_string.strip())
    return string_list
