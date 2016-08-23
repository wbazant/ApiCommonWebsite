import CheckboxTreeController from './checkboxTreeController';
import { trimBooleanQuestionAttribs } from './client/util/modelSpecificUtil';
import { getTree } from 'wdk-client/OntologyUtils';
import { isQualifying, addSearchSpecificSubtree } from 'wdk-client/CategoryUtils';
import context from './client/main';

wdk.namespace("eupathdb.attributeCheckboxTree", function(ns, $) {
  "use strict";

  // will map from summary views to attribute tree controller for that view
  let controllerMap = {};
  
  /**
   * Entry into a custom attribute selector pop-up, which appears when the user
   * clicks the Add Columns button on the header of the results table.  It has
   * been rewritten to use our React SearchableCheckboxTree component and the
   * categories ontology.
   * 
   * @param element - div from which this function was called.
   * @param attributes - attributes derived from the div - question name, record
   * class name, default selected list, current selected list, view name
   * @returns {Promise.<T>}
   */
  function setUpCheckboxTree(element, attributes) {
    let { wdkService } = context;
    let questionName = attributes.questionName;
    let recordClassName = attributes.recordClassName;
    let defaultSelectedList = attributes.defaultSelectedList.replace(/'/g,"").split(",");
    let currentSelectedList = attributes.currentSelectedList.replace(/'/g,"").split(",");
    let viewName = {'_default':'gene', 'transcript-view':'transcript'}[attributes.viewName];
    return Promise.all([
      wdkService.getOntology(),
      wdkService.findQuestion(question => question.name === questionName),
      wdkService.findRecordClass(recordClass => recordClass.name === recordClassName)
    ])
    .then(([ categoriesOntology, question, recordClass]) => {
      try {
        let categoryTree = getTree(categoriesOntology, isQualifying({
          targetType: 'attribute',
          recordClassName,
          scope: 'results'
        }));
        categoryTree = trimBooleanQuestionAttribs(question, categoryTree);
        categoryTree = addSearchSpecificSubtree(question, categoryTree);
        let selectedList = currentSelectedList || defaultSelectedList;
        let controller = new CheckboxTreeController(element, "selectedFields",
            categoryTree, selectedList, defaultSelectedList);
        controllerMap[attributes.viewName] = controller;
      }
      catch (error) {
        console.error(error);
        throw error;
      }
    })
    .catch(function(error) {
      console.error(error);
      throw new Error(error.message);
    });
  }

  function mountCheckboxTree(viewName) {
    let controller = controllerMap[viewName];
    if (controller) {
      controller.renderTree();
    }
    
  }

  function unmountCheckboxTree(viewName) {
    let controller = controllerMap[viewName];
    if (controller) {
      controller.unmountTree();
    }
  }
  
  ns.setUpCheckboxTree = setUpCheckboxTree;
  ns.mountCheckboxTree = mountCheckboxTree;
  ns.unmountCheckboxTree = unmountCheckboxTree;

});