#!/usr/bin/env python3
"""
Query PubMed to retrieve "cited by" records for a list of PubMed IDs.
Uses NCBI E-utilities API to find papers that cite given PMIDs.

The scripts makes these calls:
- Find citations with (for example, PMID 24618469):
  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pubmed&id=24618469&linkname=pubmed_pubmed_citedin&email=your.email@example.com&tool=pubmed_citation_query
- Get paper details with (if citing PMIDs were 40727939,37200190):
  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=40727939,37200190&retmode=xml&email=your.email@example.com&tool=pubmed_citation_query

Usage Examples:
  # Query citations for a specific PMID
  python query_citations.py 24618469 --output-file citations_24618469.tsv --email your.email@example.com
  
  # Query citations from a file with PMIDs
  query_citations.py --input-file immcantation-publications.tsv --output-file immcantation_citations.tsv --email your.email@example.com

  # Save corresponding author emails in the output (disabled by default)
  query_citations.py --input-file immcantation-publications.tsv --output-file immcantation_citations.tsv --email your.email@example.com --save-emails

"""

import re
import requests
import xml.etree.ElementTree as ET
import time
import csv
import argparse
import sys
import os
from datetime import datetime

class PubMedCitationQuery:
    def __init__(self, email="your.email@example.com", tool="pubmed_citation_query", save_emails=False):
        """
        Initialize the PubMed citation query tool.
        
        Args:
            email (str): Your email address (required by NCBI)
            tool (str): Name of your tool/script
            save_emails (bool): Whether to save corresponding author emails in output (default: False)
        """
        self.email = email
        self.tool = tool
        self.save_emails = save_emails
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.delay = 0.34  # NCBI recommends max 3 requests per second
        
    def search_citations(self, pmid):
        """
        Search for papers that cite the given PMID.
        
        Args:
            pmid (str): PubMed ID to search citations for
            
        Returns:
            list: List of PMIDs that cite the input PMID
        """
        # Use ELink to find citations via PMC citation database
        return self._get_citations_via_elink(pmid)
    
    def _get_citations_via_elink(self, pmid):
        """
        Use ELink to find citations via PMC citation database.
        
        Args:
            pmid (str): PubMed ID
            
        Returns:
            list: List of citing PMIDs
        """
        elink_url = f"{self.base_url}elink.fcgi"
        elink_params = {
            'dbfrom': 'pubmed',
            'db': 'pubmed', 
            'id': pmid,
            'linkname': 'pubmed_pubmed_citedin',  # Papers that cite this one
            'email': self.email,
            'tool': self.tool
        }
        
        try:
            time.sleep(self.delay)
            response = requests.get(elink_url, params=elink_params)
            
            # Check for 500 server errors specifically
            if response.status_code == 500:
                print(f"Error: NCBI server error (500) for PMID {pmid}")
                print("NCBI E-utilities is experiencing server issues. Please try again later.")
                sys.exit(1)
            
            response.raise_for_status()
            
            # Parse XML response
            root = ET.fromstring(response.content)
            citing_pmids = []
            
            # Look for linked IDs in the XML
            for link_set in root.findall('.//LinkSet'):
                for link_set_db in link_set.findall('.//LinkSetDb'):
                    for link in link_set_db.findall('.//Link'):
                        citing_id = link.find('Id')
                        if citing_id is not None:
                            citing_pmids.append(citing_id.text)
            
            return citing_pmids
            
        except requests.RequestException as e:
            # Check if it's a 500 error that wasn't caught above
            if hasattr(e, 'response') and e.response is not None and e.response.status_code == 500:
                print(f"Error: NCBI server error (500) for PMID {pmid}")
                print("NCBI E-utilities is experiencing server issues. Please try again later.")
                sys.exit(1)
            print(f"Error fetching citations for PMID {pmid}: {e}")
            return []
        except ET.ParseError as e:
            print(f"Error parsing XML response for PMID {pmid}: {e}")
            return []
    
    def get_paper_details(self, pmids):
        """
        Get detailed information for a list of PMIDs.
        
        Args:
            pmids (list): List of PubMed IDs
            
        Returns:
            list: List of dictionaries with paper details
        """
        if not pmids:
            return []
        
        all_papers = []
        batch_size = 200
        
        # Process PMIDs in batches of 200
        for i in range(0, len(pmids), batch_size):
            batch_pmids = pmids[i:i + batch_size]
            
            # EFetch to get detailed records
            efetch_url = f"{self.base_url}efetch.fcgi"
            pmid_string = ','.join(batch_pmids)
            
            efetch_params = {
                'db': 'pubmed',
                'id': pmid_string,
                'retmode': 'xml',
                'email': self.email,
                'tool': self.tool
            }
            
            try:
                time.sleep(self.delay)
                response = requests.get(efetch_url, params=efetch_params)
                
                # Check for 500 server errors specifically
                if response.status_code == 500:
                    print("Error: NCBI server error (500) while fetching paper details")
                    print("NCBI E-utilities is experiencing server issues. Please try again later.")
                    sys.exit(1)
                
                response.raise_for_status()
                
                batch_papers = self._parse_pubmed_xml(response.content)
                all_papers.extend(batch_papers)
                
                if len(pmids) > batch_size:
                    print(f"    Processed batch {i//batch_size + 1}/{(len(pmids)-1)//batch_size + 1} ({len(batch_papers)} papers)")
                
            except requests.RequestException as e:
                # Check if it's a 500 error that wasn't caught above
                if hasattr(e, 'response') and e.response is not None and e.response.status_code == 500:
                    print("Error: NCBI server error (500) while fetching paper details")
                    print("NCBI E-utilities is experiencing server issues. Please try again later.")
                    sys.exit(1)
                print(f"Error fetching details for batch {i//batch_size + 1}: {e}")
                continue
        
        return all_papers
    
    def _parse_pubmed_xml(self, xml_content):
        """
        Parse PubMed XML response to extract paper details.
        
        Args:
            xml_content (bytes): XML content from PubMed
            
        Returns:
            list: List of paper detail dictionaries
        """
        papers = []
        
        try:
            root = ET.fromstring(xml_content)
            
            for article in root.findall('.//PubmedArticle'):
                paper = {
                    'pmid': '',
                    'title': '',
                    'authors': '',
                    'corresponding_authors': '',
                    'corresponding_emails': '',
                    'journal': '',
                    'year': '',
                    'publication_date': '',
                    'abstract': '',
                    'doi': '',
                    'pmcid': ''
                }
                
                # Get PMID
                pmid_elem = article.find('.//PMID')
                if pmid_elem is not None:
                    paper['pmid'] = pmid_elem.text
                
                # Get title
                title_elem = article.find('.//ArticleTitle')
                if title_elem is not None:
                    paper['title'] = ''.join(title_elem.itertext()).strip()
                
                # Get authors and corresponding authors
                authors = []
                corresponding_authors = []
                corresponding_emails = []
                all_emails_in_article = []
                
                # First pass: collect all authors
                author_list = []
                for author in article.findall('.//Author'):
                    last_name = author.find('LastName')
                    first_name = author.find('ForeName')
                    if last_name is not None:
                        author_name = last_name.text
                        if first_name is not None:
                            author_name = f"{last_name.text} {first_name.text}"
                        
                        author_info = {
                            'name': author_name,
                            'element': author,
                            'emails': [],
                            'is_corresponding': False
                        }
                        
                        # Check for emails in affiliation
                        for affiliation_info in author.findall('.//AffiliationInfo'):
                            affiliation_text = affiliation_info.find('Affiliation')
                            if affiliation_text is not None and affiliation_text.text:
                                email_pattern = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b'
                                email_matches = re.findall(email_pattern, affiliation_text.text)
                                author_info['emails'].extend(email_matches)
                                all_emails_in_article.extend(email_matches)
                        
                        # Check if explicitly marked as corresponding author
                        # Look for attributes or special elements that indicate corresponding author
                        author_attr = author.get('ValidYN')
                        if author_attr == 'Y':
                            # Some records mark corresponding authors this way
                            author_info['is_corresponding'] = True
                        
                        author_list.append(author_info)
                        authors.append(author_name)
                
                # Second pass: determine corresponding authors
                # Strategy 1: Authors with emails are corresponding authors
                authors_with_emails = [auth for auth in author_list if auth['emails']]
                
                if authors_with_emails:
                    # Use authors who have emails as corresponding authors
                    for auth in authors_with_emails:
                        unique_emails = list(dict.fromkeys(auth['emails']))
                        email_str = ', '.join(unique_emails)
                        corresponding_authors.append(f"{auth['name']} <{email_str}>")
                        corresponding_emails.extend(unique_emails)
                
                else:
                    # Strategy 2: Check for emails mentioned elsewhere in the article
                    # Look for contact information or author notes
                    contact_info = []
                    
                    # Check for email in other parts of the article
                    for elem in article.findall('.//*'):
                        if elem.text and '@' in elem.text:
                            email_pattern = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b'
                            emails_in_text = re.findall(email_pattern, elem.text)
                            contact_info.extend(emails_in_text)
                    
                    if contact_info and author_list:
                        # If we found emails but couldn't associate them with specific authors,
                        # assume the last author is corresponding (common convention)
                        last_author = author_list[-1]
                        unique_emails = list(dict.fromkeys(contact_info))
                        email_str = ', '.join(unique_emails)
                        corresponding_authors.append(f"{last_author['name']} <{email_str}>")
                        corresponding_emails.extend(unique_emails)
                    
                    elif author_list:
                        # Strategy 3: If no emails found anywhere, assume last author is corresponding
                        # (common academic convention) but don't include email
                        last_author = author_list[-1]
                        corresponding_authors.append(last_author['name'])
                
                # Remove duplicate emails while preserving order
                corresponding_emails = list(dict.fromkeys(corresponding_emails))
                
                paper['authors'] = ', '.join(authors)
                paper['corresponding_authors'] = '; '.join(corresponding_authors)
                paper['corresponding_emails'] = '; '.join(list(dict.fromkeys(corresponding_emails)))
                
                # Get abstract
                abstract_parts = []
                for abstract in article.findall('.//Abstract'):
                    for abstract_text in abstract.findall('.//AbstractText'):
                        label = abstract_text.get('Label', '')
                        text = ''.join(abstract_text.itertext()).strip()
                        if label:
                            abstract_parts.append(f"{label}: {text}")
                        else:
                            abstract_parts.append(text)
                
                paper['abstract'] = ' '.join(abstract_parts)
                
                # Get journal
                journal_elem = article.find('.//Journal/Title')
                if journal_elem is not None:
                    paper['journal'] = journal_elem.text
                
                # Get year and publication date
                year_elem = article.find('.//PubDate/Year')
                month_elem = article.find('.//PubDate/Month')
                day_elem = article.find('.//PubDate/Day')
                
                if year_elem is not None:
                    paper['year'] = year_elem.text
                    
                    # Format publication date as YYYY-MM-DD
                    year = year_elem.text
                    month = '01'  # Default to January if no month
                    day = '01'    # Default to 1st if no day
                    
                    if month_elem is not None:
                        month_text = month_elem.text
                        # Handle month names (convert to numbers)
                        month_map = {
                            'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04',
                            'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08',
                            'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'
                        }
                        if month_text in month_map:
                            month = month_map[month_text]
                        elif month_text.isdigit() and 1 <= int(month_text) <= 12:
                            month = f"{int(month_text):02d}"
                    
                    if day_elem is not None and day_elem.text.isdigit():
                        day = f"{int(day_elem.text):02d}"
                    
                    paper['publication_date'] = f"{year}-{month}-{day}"
                
                # Get DOI
                for article_id in article.findall('.//ArticleId'):
                    if article_id.get('IdType') == 'doi':
                        paper['doi'] = article_id.text
                        break
                
                # Get PMCID
                for article_id in article.findall('.//ArticleId'):
                    if article_id.get('IdType') == 'pmc' and article_id.text:
                        pmcid_value = article_id.text.strip()
                        # Only accept PMCIDs that start with "PMC" (valid format)
                        if pmcid_value.startswith('PMC'):
                            paper['pmcid'] = pmcid_value
                            break
                
                papers.append(paper)
                
        except ET.ParseError as e:
            print(f"Error parsing XML: {e}")
            
        return papers
    
    def process_pmid_list(self, pmids, output_file=None, titles=None):
        """
        Process a list of PMIDs and get their citations.
        
        Args:
            pmids (list): List of PubMed IDs to process
            output_file (str): Optional output file path
            titles (dict): Optional dict mapping PMID to title for display
            
        Returns:
            dict: Dictionary mapping input PMIDs to citing papers
        """
        results = {}
        
        for i, pmid in enumerate(pmids):
            title_preview = ""
            if titles and pmid in titles:
                t = titles[pmid]
                title_preview = f" - {t[:60]}..." if len(t) > 60 else f" - {t}"
            print(f"Processing PMID {pmid} ({i+1}/{len(pmids)}){title_preview}...")
            
            # Get citing PMIDs
            citing_pmids = self.search_citations(pmid)
            print(f"  Found {len(citing_pmids)} citing papers")
            
            if citing_pmids:
                # Get details for citing papers
                citing_papers = self.get_paper_details(citing_pmids)
                if len(citing_papers) < len(citing_pmids):
                    print(f"  WARNING: Retrieved details for only {len(citing_papers)}/{len(citing_pmids)} citing papers (possible fetch failure)")
                results[pmid] = citing_papers
            else:
                results[pmid] = []
            
            # Be nice to NCBI servers
            time.sleep(self.delay)
        
        if output_file:
            self._save_results(results, output_file)
            
        return results
    
    def _save_results(self, results, output_file):
        """
        Save results to a TSV file.
        
        Args:
            results (dict): Results dictionary
            output_file (str): Output file path
        """
        fieldnames = ['original_pmid', 'citing_pmid', 'title', 'authors', 'corresponding_authors', 
                     'corresponding_emails', 'journal', 'year', 'publication_date', 'abstract', 'doi', 'pmcid']

        def _strip_emails(text):
            """Remove <email> portions from a corresponding-author string."""
            return re.sub(r'\s*<[^>]+>', '', text).strip()
        
        # Get current date for the download timestamp
        download_date = datetime.now().strftime("%Y-%m-%d")
        
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            # Write download date as a comment line at the top
            f.write(f"# Downloaded: {download_date}\n")
            
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', quoting=csv.QUOTE_ALL)
            writer.writeheader()
            
            for original_pmid, citing_papers in results.items():
                if not citing_papers:
                    # Write a row even if no citations found
                    writer.writerow({
                        'original_pmid': original_pmid,
                        'citing_pmid': '',
                        'title': '',
                        'authors': '',
                        'corresponding_authors': '',
                        'corresponding_emails': '',
                        'journal': '',
                        'year': '',
                        'publication_date': '',
                        'abstract': '',
                        'doi': '',
                        'pmcid': ''
                    })
                else:
                    for paper in citing_papers:
                        corresponding_authors = paper['corresponding_authors']
                        corresponding_emails = paper['corresponding_emails']
                        if not self.save_emails:
                            corresponding_authors = _strip_emails(corresponding_authors)
                            corresponding_emails = ''
                        writer.writerow({
                            'original_pmid': original_pmid,
                            'citing_pmid': paper['pmid'],
                            'title': paper['title'],
                            'authors': paper['authors'],
                            'corresponding_authors': corresponding_authors,
                            'corresponding_emails': corresponding_emails,
                            'journal': paper['journal'],
                            'year': paper['year'],
                            'publication_date': paper['publication_date'],
                            'abstract': paper['abstract'],
                            'doi': paper['doi'],
                            'pmcid': paper['pmcid']
                        })

def main():
    parser = argparse.ArgumentParser(description='Query PubMed for citation data')
    parser.add_argument('pmids', nargs='*', help='PubMed IDs to query (space-separated)')
    parser.add_argument('--input-file', '-i', 
                       help='TSV/CSV file with a column named "pmid" or "PMID"')
    parser.add_argument('--output-file', '-o', help='Output TSV file')
    parser.add_argument('--email', '-e', default='your.email@example.com', 
                       help='Your email address (required by NCBI)')
    parser.add_argument('--save-emails', action='store_true', default=False,
                       help='Save corresponding author emails in the output TSV (default: False)')
    
    args = parser.parse_args()
    
    # Collect PMIDs from command line arguments and/or input file
    pmids = []
    titles = {}  # optional mapping of pmid -> title
    
    if args.pmids:
        pmids.extend(args.pmids)
    
    if args.input_file:
        try:
            # Try to read as TSV/CSV with headers
            file_ext = os.path.splitext(args.input_file)[1].lower()
            delimiter = '\t' if file_ext == '.tsv' else ','
            
            with open(args.input_file, 'r', encoding='utf-8') as f:
                # Read the entire content first
                content = f.read()
                f.seek(0)
                
                # First, try to read as a structured file with headers
                try:
                    reader = csv.DictReader(f, delimiter=delimiter)
                    headers = reader.fieldnames
                    
                    # Look for PMID column (case insensitive)
                    pmid_column = None
                    if headers:  # Check if headers exist
                        for header in headers:
                            if header.lower() in ['pmid', 'pubmed_id', 'pubmedid']:
                                pmid_column = header
                                break
                    
                    if pmid_column:
                        # Look for a title column (case insensitive)
                        title_column = None
                        if headers:
                            for header in headers:
                                if header.lower() == 'title':
                                    title_column = header
                                    break
                        # Read PMIDs (and titles if available) from the structured file
                        for row in reader:
                            pmid = row.get(pmid_column, '').strip()
                            if pmid and pmid.isdigit():
                                pmids.append(pmid)
                                if title_column:
                                    titles[pmid] = row.get(title_column, '').strip()
                        print(f"Loaded {len([p for p in pmids if p])} PMIDs from column '{pmid_column}' in {args.input_file}")
                    else:
                        # Fallback: treat as plain text file with one PMID per line
                        for line in content.splitlines():
                            pmid = line.strip()
                            if pmid and pmid.isdigit():
                                pmids.append(pmid)
                        print(f"Loaded {len([p for p in pmids if p])} PMIDs from {args.input_file} (plain text format)")
                        
                except csv.Error:
                    # If CSV parsing fails, treat as plain text
                    for line in content.splitlines():
                        pmid = line.strip()
                        if pmid and pmid.isdigit():
                            pmids.append(pmid)
                    print(f"Loaded {len([p for p in pmids if p])} PMIDs from {args.input_file} (plain text format)")
                            
        except FileNotFoundError:
            print(f"Error: Input file {args.input_file} not found")
            sys.exit(1)
        except Exception as e:
            print(f"Error reading {args.input_file}: {e}")
            sys.exit(1)
    
    if not pmids:
        print("Error: No PMIDs provided. Use --help for usage information.")
        sys.exit(1)
    
    # Check if output directory exists when output file is specified
    if args.output_file:
        output_dir = os.path.dirname(args.output_file)
        if output_dir and not os.path.exists(output_dir):
            print(f"Error: Output directory '{output_dir}' does not exist")
            print("Please create the directory first or specify a different output path")
            sys.exit(1)
    
    # Remove duplicates while preserving order
    pmids = list(dict.fromkeys(pmids))
    
    print(f"Processing {len(pmids)} PMIDs...")
    
    # Initialize the query tool
    query_tool = PubMedCitationQuery(email=args.email, save_emails=args.save_emails)
    
    # Process the PMIDs
    results = query_tool.process_pmid_list(pmids, args.output_file, titles=titles)
    
    # Print summary
    total_citations = sum(len(citing_papers) for citing_papers in results.values())
    unique_citations = len({p['pmid'] for papers in results.values() for p in papers})
    print("\nSummary:")
    print(f"  Processed {len(pmids)} PMIDs")
    print(f"  Total citations: {total_citations} ({unique_citations} unique citations)")
    print("")
    print(f"  {'PMID':<12} {'Citing papers':<15} Title")
    print(f"  {'-'*12} {'-'*15} {'-'*40}")
    for original_pmid, citing_papers in results.items():
        title_str = titles.get(original_pmid, "")
        if len(title_str) > 60:
            title_str = title_str[:60] + "..."
        print(f"  {original_pmid:<12} {len(citing_papers):<15} {title_str}")

    if args.output_file:
        print(f"\n  Results saved to {args.output_file}")

if __name__ == '__main__':
    main()
