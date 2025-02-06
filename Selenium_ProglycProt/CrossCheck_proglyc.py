import sys

import numpy as np
import selenium
from selenium import webdriver
from selenium import webdriver
from selenium.webdriver import ActionChains
from selenium.webdriver.chrome.service import Service as ChromeService
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait
from webdriver_manager.chrome import ChromeDriverManager
import time
from selenium.webdriver.support.ui import Select
import pandas as pd


driver = webdriver.Chrome(service=ChromeService(ChromeDriverManager().install()))
wait = WebDriverWait(driver, 10)
actions = ActionChains(driver)


def check_if(value1, checks):
    driver.get("https://www.proglycprot.org")

    cssmenu = driver.find_element(by=By.ID, value='cssmenu')
    nav_list = cssmenu.find_element(by=By.TAG_NAME, value='ul')
    nav_second_element = nav_list.find_elements(by=By.TAG_NAME, value='li')[1]
    actions.move_to_element(nav_second_element).perform()

    wait.until(EC.visibility_of_element_located((By.XPATH, "//a[@data-target='#search']"))).click()
    m = wait.until(EC.visibility_of_element_located((By.ID, 'login-form')))

    Select(m.find_element(by=By.ID, value='searchcriteria1')).select_by_value('uniprotid')
    time.sleep(1)
    try:
        Select(m.find_element(by=By.ID, value='value1')).select_by_visible_text(value1)
    except:
        return "NOT IN DATABASE"
    m.submit()

    out = wait.until(EC.presence_of_element_located((By.ID, "columns")))
    table = out.find_element(by=By.ID, value="no-more-tables")

    table.find_element(by=By.TAG_NAME, value='a').click()

    tt = wait.until(EC.visibility_of_element_located((By.XPATH, "//table[@id='customers']")))

    all_rows = tt.find_elements(By.TAG_NAME, "tr")
    for row in all_rows:
        td = row.find_elements(By.TAG_NAME, "td")
        if td[0].text == "UniProtKB Sequence":
            href = td[1].find_element(By.TAG_NAME, "a").get_attribute("href")
            driver.get(href)
            break

    try:
        sequence = "".join(driver.find_element(by=By.TAG_NAME, value="pre").text.split('\n')[1:])
    except:
        return "SEQUENCE IS EMPTY"

    pos = {int(check[1:]): check[0] for check in checks.split(", ")}
    return str({f"{pos[i]}{i}": sequence[i-1] == pos[i] for i in pos.keys()})


seq_empty = []
not_in_db = []
df = pd.read_csv('check_output_file.csv')
for ii in df.index:
    row = df.loc[ii]
    if not pd.isna(row['check_output']):
        continue
    try:
        returnvalue = check_if(value1=row['unit_pro_kb'], checks=row['glycosites'])
    except:
        print('Something wrong with: ', row['unit_pro_kb'])
        continue
    if returnvalue == 'SEQUENCE IS EMPTY':
        seq_empty.append(row['unit_pro_kb'])
        continue
    elif returnvalue == "NOT IN DATABASE":
        not_in_db.append(row['unit_pro_kb'])
        continue

    df.loc[ii, 'check_output'] = returnvalue
    df.to_csv('check_output_file.csv', index=False)

driver.close()

print(not_in_db)
print(seq_empty)

not_in_db = ['WP065336528', 'A1VZX2', 'A0A0P1D0G7', 'YP1778551', 'Q9AET1', 'Q2N2Q5', 'P9WIR6', 'YP1778491', 'YP1777381', 'UPI00062B1A4A', 'Q46080', 'C0X1N7', 'A0A0U0B015', 'NP2184071', 'NP2182801', 'NP2174571', 'NP2173151', 'NP2172191', 'NP2168061', 'NP2166611', 'NP2165431', 'NP2161931', 'NP2161571', 'NP2160431', 'NP2157691', 'NP2154711', 'NP2152201', 'NP2149361', 'NP2148641', 'NP2148621', 'NP2148251', 'NP2147301', 'NP2145801', 'NP2145341', 'CCP460921', 'CCP460671', 'CCP460461', 'CCP450211', 'CCP449981', 'CCP441431', 'CCP438521', 'CCP438481', 'CCP432731', 'CAD554621', 'CAD309531', 'CAC443221', 'CAB930671', 'CAB884741', 'CAB773501', 'CAB655681', 'CAB408571', 'CAB408531', 'CAB384761', 'CAA223621', 'CAA198521', 'CAA185161', 'CAA158821', 'ANZ820431', 'ANZ819371', 'ABO116231', 'AAF134001']
seq_empty = ['C3GCL2', 'Q5F718', 'A0R2Q4', 'Q0PAL6']
