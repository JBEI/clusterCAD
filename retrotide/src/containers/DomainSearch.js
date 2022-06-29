import React from 'react';
import { connect } from 'react-redux';
import { beginDomainSearch } from '../redux/actions/actions';
import {clusterCADDomainSearch} from '../redux/middleware/api';
import Button from '../components/Button';
import ModuleBuilder from '../components/ModuleBuilder'

const mapDispatchToProps = dispatch => {
  return {
    beginDomainSearch: modules => dispatch(beginDomainSearch(modules)),
    dispatch,
  }
};

class DomainSearch extends React.Component {

  constructor(props) {
    super(props);

    const LoadingPKSList = {
      KS:  {domainName: 'KS', present: true},
      AT:  {domainName: 'AT', present: true, options: {
        substrate: ['mal', 'methylmal', 'trans'],
        selected: 'mal',
      }},
      ACP: {domainName: 'ACP', present: true},       
    };

    const TerminatingPKSList = {
      KS:  {domainName: 'KS', present: true},
      AT:  {domainName: 'AT', present: true, options: {
        substrate: ['mal', 'methylmal', 'trans'],
        selected: 'mal',
      }},
      DH:  {domainName: 'DH', present: false},
      ER:  {domainName: 'ER', present: false},
      KR:  {domainName: 'KR', present: false, options: {
        stereochemistry: ['A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'any'],
        selected: 'any',
      }},
      ACP: {domainName: 'ACP', present: true}, 
      TE:  {domainName: 'TE', present: true},        
    };

    const ButtonList = {
      KR: {domainName: 'KR', domains: ['KR'], disabled: false},
      DH_KR: {domainName: 'DH-KR', domains: ['DH', 'KR'], disabled: false},
      DH_ER_KR: {domainName: 'DH-ER-KR', domains: ['DH', 'ER', 'KR'], disabled: false},
    }

    this.state = {
      ModuleArray: [],
      LoadingModule: {
        key: 'loading-0',
        type: 'loading',
        domainList: LoadingPKSList,
      },
      TerminatingModule: {
        key: 'terminating-n',
        type: 'terminating',
        domainList: TerminatingPKSList,
        buttonList: ButtonList,
      },
      DomainList: 'PKS',
    }
  }

  buildModules = (newLength) => {
    let PKSDomainList = {
      KS:  {domainName: 'KS', present: true},
      AT:  {domainName: 'AT', present: true, options: {
        substrate: ['mal', 'methylmal', 'trans'],
        selected: 'mal',
      }},
      DH:  {domainName: 'DH', present: false},
      ER:  {domainName: 'ER', present: false},
      KR:  {domainName: 'KR', present: false, options: {
        stereochemistry: ['A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'any'],
        selected: 'any',
      }},
      ACP: {domainName: 'ACP', present: true}, 
    };
    let ButtonList = {
      KR: {domainName: 'KR', domains: ['KR'], disabled: false},
      DH_KR: {domainName: 'DH-KR', domains: ['DH', 'KR'], disabled: false},
      DH_ER_KR: {domainName: 'DH-ER-KR', domains: ['DH', 'ER', 'KR'], disabled: false},
    }
    let defaultArray = Array.from(
      { length: newLength },
      (x, index) => ({
        index: index,
        key: 'extending-' + index,
        domainList: PKSDomainList,
        buttonList: ButtonList,
        type: 'extending',
      })
    );
    this.setState({ ModuleArray: defaultArray });
  }

  addModule = () => {
    let currentPlusOne = this.state.ModuleArray.length + 1;
    this.buildModules(currentPlusOne);
  }

  parseModuleObject = (module, index) => {
    let {key, type, domainList, buttonList} = module;
    console.log(index);
    return (
      <ModuleBuilder 
        index={index} 
        key={key}
        id={key}
        deleteFunction={this.deleteModule}
        updateFunction={this.updateModule}
        domainList={domainList}
        buttonList={buttonList}
        type={type} />
    );
  }

  deleteModule = (moduleKey) => {
    let currentModules = this.state.ModuleArray;
    let moduleIndex = currentModules.findIndex((modObj) => modObj.key === moduleKey);
    if(moduleIndex >= 0) {
      currentModules.splice(moduleIndex, 1);
    }
    this.setState({
      ModuleArray: currentModules,
    });
  }

  updateModule = (moduleKey, newModuleContent) => {
    let currentModules = this.state.ModuleArray;
    let moduleIndex = currentModules.findIndex((modObj) => modObj.key === moduleKey);
    let moduleToUpdate = currentModules[moduleIndex];
    currentModules[moduleIndex] = {... moduleToUpdate, DomainList: newModuleContent};
    this.setState({
      ModuleArray: currentModules
    });
  }

  submitSearch = () => {
    let loading = this.state.LoadingModule;
    let terminating = this.state.TerminatingModule;
    let currentModules = this.state.ModuleArray;
    // push and unshift are mutators
    // this works fine even if ModuleArray has length 0
    currentModules.unshift(loading);
    currentModules.push(terminating);
    // dispatch action
    this.props.beginDomainSearch(currentModules);
    // call async function. is this where this goes?
    clusterCADDomainSearch(currentModules);
  }

  render() {
    const ExtendingArray = this.state.ModuleArray;
    return (
      <div className='DomainSearch form'>
        <h3>Construct Modules</h3>
        <Button onClick={() => { this.addModule() }}> Add Module + </Button>
        <Button onClick={() => { this.submitSearch() }} className="submit"> Submit </Button>
        <div className="ModuleListWrapper">
          { this.parseModuleObject(this.state.LoadingModule, -1) }
          { (ExtendingArray && ExtendingArray.length > 0) &&
            ExtendingArray.map((exModule, index) => {
              return (
                this.parseModuleObject(exModule, index)
              )
            }) 
          }
          { this.parseModuleObject(this.state.TerminatingModule, ExtendingArray ? ExtendingArray.length : 1) }
        </div>
      </div>
    )
  }
};

export default connect(null, mapDispatchToProps)(DomainSearch);