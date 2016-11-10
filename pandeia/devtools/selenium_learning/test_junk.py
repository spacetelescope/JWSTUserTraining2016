# -*- coding: utf-8 -*-
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.common.exceptions import NoSuchElementException
from selenium.common.exceptions import NoAlertPresentException
import unittest, time, re

class SeleniumIdeTestCases(unittest.TestCase):
    def setUp(self):
        self.driver = webdriver.Firefox()
        self.driver.implicitly_wait(30)
        self.base_url = "http://etcbrady.stsci.edu:4242/"
        self.verificationErrors = []
        self.accept_next_alert = True
    
    def test_selenium_ide_test_cases(self):
        driver = self.driver
        driver.find_element_by_link_text("Instrument Setup").click()
        driver.find_element_by_xpath("//button[@onclick='recompute();']").click()
        driver.find_element_by_xpath("//tr[@id='1']/td[4]").click()
        driver.find_element_by_link_text("Focal Plane Rate").click()
        driver.find_element_by_link_text("Target").click()
        driver.find_element_by_link_text("Detector").click()
        driver.find_element_by_xpath("(//input[@id='overplotcheck'])[5]").click()
        driver.find_element_by_link_text("2D SNR").click()
        driver.find_element_by_link_text("Detector").click()
        driver.find_element_by_link_text("Focal Plane Rate").click()
        driver.find_element_by_link_text("Report").click()
        driver.find_element_by_xpath("//tr[@id='5']/td[4]").click()
        driver.find_element_by_xpath("(//input[@id='overplotcheck'])[5]").click()
        driver.find_element_by_id("overplotcheck").click()
        driver.find_element_by_link_text("Detector Setup").click()
        driver.find_element_by_xpath("//button[@onclick='recompute();']").click()
    
    def is_element_present(self, how, what):
        try: self.driver.find_element(by=how, value=what)
        except NoSuchElementException, e: return False
        return True
    
    def is_alert_present(self):
        try: self.driver.switch_to_alert()
        except NoAlertPresentException, e: return False
        return True
    
    def close_alert_and_get_its_text(self):
        try:
            alert = self.driver.switch_to_alert()
            alert_text = alert.text
            if self.accept_next_alert:
                alert.accept()
            else:
                alert.dismiss()
            return alert_text
        finally: self.accept_next_alert = True
    
    def tearDown(self):
        self.driver.quit()
        self.assertEqual([], self.verificationErrors)

if __name__ == "__main__":
    unittest.main()
