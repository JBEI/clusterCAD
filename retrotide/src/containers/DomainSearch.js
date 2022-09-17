import React from 'react';
import { connect } from 'react-redux';
import { beginDomainSearch } from '../redux/actions/actions';
import {clusterCADDomainSearch} from '../redux/middleware/api';
import { Redirect } from 'react-router-dom';
import Button from '../components/Button';
import ModuleBuilder from '../components/ModuleBuilder';
import addIcon from '../images/add-circle-fill.png';

const mapDispatchToProps = dispatch => {
  return {
    clusterCADDomainSearch: modules => dispatch(clusterCADDomainSearch(modules, token)),
    dispatch,
  }
};

const cookie = document.cookie;
const csrftoken = cookie.match(/csrftoken=([a-zA-Z0-9]+);?/);

class DomainSearch extends React.Component {

  constructor(props) {
    super(props);

    // all of these consts need to go into a json file or something
    // they should come in via a reducer that will set them up in this
    // format so we can use them. We need to be able to parse either PKS
    // or potentially other formats

    // right now these aren't being stored in state properly

    const LoadingPKSList = {
      KS:  {domainName: 'KS', present: true},
      AT:  {domainName: 'AT', present: true, options: {
        substrate: ['mal', 'mmal', 'mxmal', 'mxmal_ACP', 'emal', 'allylmal', 'butmal', 'hmal', 
          'isobutmal', 'D-isobutmal', 'DCP', 'hexmal', '2-oxobutmal', '3-me-hexmal'],
        selected: 'mal',
      }},
      ACP: {domainName: 'ACP', present: true},       
    };

    const TerminatingPKSList = {
      KS:  {domainName: 'KS', present: true},
      AT:  {domainName: 'AT', present: true, options: {
        substrate: ['mal', 'mmal', 'mxmal', 'mxmal_ACP', 'emal', 'allylmal', 'butmal', 'hmal', 
          'isobutmal', 'D-isobutmal', 'DCP', 'hexmal', '2-oxobutmal', '3-me-hexmal'],
        selected: 'mal',
      }},
      DH:  {domainName: 'DH', present: false},
      ER:  {domainName: 'ER', present: false},
      KR:  {domainName: 'KR', present: false, options: {
        stereochemistry: ['A1', 'A2', 'B1', 'B2', 'C1', 'C2'],
        selected: 'B1',
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
        DomainList: LoadingPKSList,
      },
      TerminatingModule: {
        key: 'terminating-n',
        type: 'terminating',
        DomainList: TerminatingPKSList,
        ButtonList: ButtonList,
      },
      DomainList: 'PKS',
      redirectToResults: false,
    }
  }

  // same thing here, these lists need to be moved out and brought into local state properly

  buildModules = (newLength) => {
    let PKSDomainList = {
      KS:  {domainName: 'KS', present: true},
      AT:  {domainName: 'AT', present: true, options: {
        substrate: ['mal', 'mmal', 'mxmal', 'mxmal_ACP', 'emal', 'allylmal', 'butmal', 'hmal', 
          'isobutmal', 'D-isobutmal', 'DCP', 'hexmal', '2-oxobutmal', '3-me-hexmal'],
        selected: 'mal',
      }},
      DH:  {domainName: 'DH', present: false},
      ER:  {domainName: 'ER', present: false},
      KR:  {domainName: 'KR', present: false, options: {
        stereochemistry: ['A1', 'A2', 'B1', 'B2', 'C1', 'C2'],
        selected: 'B1',
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
        DomainList: PKSDomainList,
        ButtonList: ButtonList,
        type: 'extending',
      })
    );
    this.setState({ ModuleArray: defaultArray });
  }

  // runs when you hit the add module button
  addModule = () => {
    let currentPlusOne = this.state.ModuleArray.length + 1;
    this.buildModules(currentPlusOne);
  }

  // this takes modules from the list and returns the JSX we need to put them in the dom
  parseModuleObject = (module, index) => {
    let {key, type, DomainList, ButtonList} = module;
    return (
      <ModuleBuilder 
        index={index} 
        key={key}
        id={key}
        deleteFunction={this.deleteModule}
        updateFunction={this.updateModule}
        DomainList={DomainList}
        ButtonList={ButtonList}
        type={type} />
    );
  }

  // this function is passed down to the module class
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

  // this function is passed down to the module class
  // the container needs to know when updates hapen because it holds the master list
  // of modules to submit to the backend
  updateModule = (moduleKey, newModuleContent) => {
    let loading = this.state.LoadingModule;
    let terminating = this.state.TerminatingModule;
    let currentModules = this.state.ModuleArray;
    // extending loading and terminating are stored separately so we check each one
    if(moduleKey == 'loading-0') {
      let newLoading = {... loading, DomainList: newModuleContent};
      this.setState({
        LoadingModule: newLoading,
      });
    } else if(moduleKey == 'terminating-n') {
      let newTerminating = {... terminating, DomainList: newModuleContent};
      this.setState({
        TerminatingModule: newTerminating,
      });
    } else {
      let moduleIndex = currentModules.findIndex((modObj) => modObj.key === moduleKey);
      if(moduleIndex > -1) {
        let moduleToUpdate = currentModules[moduleIndex];
        currentModules[moduleIndex] = {... moduleToUpdate, DomainList: newModuleContent};
        this.setState({
          ModuleArray: currentModules,
        });
      } else {
        console.log("Module to update not found! " + moduleKey);
      }
    }
  }

  // the async call to the backend. This is currently communicating with the backend
  submitSearch = () => {
    let loading = this.state.LoadingModule;
    let terminating = this.state.TerminatingModule;
    let currentModules = this.state.ModuleArray;
    let token = csrftoken[1] || "token not found"; // the group with just the token string is at index 1
    // push and unshift are mutators
    // this works fine even if ModuleArray has length 0
    currentModules.unshift(loading);
    currentModules.push(terminating);
    // dispatch actions
    beginDomainSearch(currentModules);
    // call async function. is this where this goes?
    this.props.dispatch(clusterCADDomainSearch(currentModules, token));
    setTimeout(() => {this.setState({redirectToResults: true})}, 1000);
  }

  // we need to pass in the current array so we can use its length to set the indices of modules correctly
  // maybe theres a nicer way of doing that?
  render() {
    const ExtendingArray = this.state.ModuleArray;
    return (
      <div className='DomainSearch form'>
        { this.state.redirectToResults ? <Redirect to='/searchResults' /> :
          < >
            <h3>Domain Architecture Search</h3>
            <Button onClick={() => { this.addModule() }} className="AddDomain"> Add Module <img src={addIcon} /> </Button>
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
          </ >
        }
      </div>
    )
  }
};

export default connect(null, mapDispatchToProps)(DomainSearch);
