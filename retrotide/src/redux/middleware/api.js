// async actions and calls to the client / backend go here
// we need access to the client if any
// redux-thunk is imported toplevel so we should be able to access it and dispatch from here but we'll see

const clusterCADSeqSearch = (molecule, token) => {
  console.log('hit api seqsearch with ' + molecule);
  return async function beginMoleculeSearch(dispatch, getState) {
    const requestObject = { // object structure needed for clusterCAD request, probably shouldn't be here
      // csrfmiddlewaretoken: "",
      // sdf: molecule,
      // draw: 1,
      // cutoff: 0.0,
      // maxCompounds: 10,
      integer: 0;
    };
    const response = await clusterCAD.post('/api/structureSearch/', {requestObject}); // callback
    dispatch({ type: 'jobAdded', payload: {response} });
  }
}

// export function synchWrapper (someInput) {
//   return async function beginMoleculeSearch(dispatch, getState) {
//     // example
//     const response = await client.post('/fakeApi/todos', { todo: initialTodo });
//     dispatch({ type: 'todos/todoAdded', payload: response.todo });
//   }
// }

export { clusterCADSeqSearch };
// then import this function in the component, which has no idea it's async
// remember to update reducers if dispatching an action